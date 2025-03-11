#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <cctype>
#include <stdexcept>

namespace fs = std::filesystem;

// Helper function to trim whitespace from both ends of a string.
std::string trim(const std::string &s) {
    size_t start = s.find_first_not_of(" \t\n\r");
    if (start == std::string::npos)
        return "";
    size_t end = s.find_last_not_of(" \t\n\r");
    return s.substr(start, end - start + 1);
}

// Helper function to split a string by whitespace.
std::vector<std::string> split(const std::string &line) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss >> token)
        tokens.push_back(token);
    return tokens;
}

// Joins a vector of strings into one string with a space separator.
std::string join(const std::vector<std::string> &tokens) {
    std::ostringstream oss;
    for (size_t i = 0; i < tokens.size(); i++) {
        if (i > 0)
            oss << " ";
        oss << tokens[i];
    }
    return oss.str();
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: merge_axt <directory> [output_file]" << std::endl;
        return 1;
    }

    std::string dirPath = argv[1];
    std::vector<std::string> mergedLines;

    // Iterate over all files in the provided directory (non-recursive).
    try {
        for (const auto &entry : fs::directory_iterator(dirPath)) {
            if (entry.is_regular_file()) {
                std::string ext = entry.path().extension().string();
                if (ext == ".axt") { // process only .axt files
                    std::ifstream infile(entry.path());
                    if (!infile) {
                        std::cerr << "Error opening file: " << entry.path() << std::endl;
                        continue;
                    }
                    std::string line;
                    while (std::getline(infile, line)) {
                        mergedLines.push_back(line);
                    }
                }
            }
        }
    } catch (const fs::filesystem_error &e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return 1;
    }

    // Clean mergedLines: remove lines that start with "##" unless they are among the first three lines.
    std::vector<std::string> cleanedLines;
    for (size_t i = 0; i < mergedLines.size(); i++) {
        // Always keep the first three lines
        if (i < 3) {
            cleanedLines.push_back(mergedLines[i]);
        } else {
            std::string trimmedLine = trim(mergedLines[i]);
            if (trimmedLine.size() >= 2 && trimmedLine.substr(0,2) == "##")
                continue;
            cleanedLines.push_back(mergedLines[i]);
        }
    }

    // Make sure there are at least three header lines.
    if (cleanedLines.size() < 3) {
        std::cerr << "Error: merged file has fewer than three lines of header." << std::endl;
        return 1;
    }

    // The first three lines are header.
    std::vector<std::string> header(cleanedLines.begin(), cleanedLines.begin() + 3);

    // Process the rest of the file into alignment blocks (each with 3 non-empty lines).
    std::vector< std::vector<std::string> > alignments;
    std::vector<std::string> block;
    for (size_t i = 3; i < cleanedLines.size(); i++) {
        std::string line = trim(cleanedLines[i]);
        if (line.empty())
            continue;
        block.push_back(cleanedLines[i]);
        if (block.size() == 3) {
            alignments.push_back(block);
            block.clear();
        }
    }
    if (!block.empty()) {
        std::cerr << "Warning: Incomplete alignment block found (not a multiple of 3 non-empty lines)." << std::endl;
    }

    // Sort the alignment blocks.
    // The sort key is derived from the first line of each block:
    // Compare by the second element (target chromosome),
    // then by the third element (target start as integer),
    // then by the fifth element (query chromosome),
    // and finally by the sixth element (query start as integer).
    std::sort(alignments.begin(), alignments.end(), [](const std::vector<std::string>& a, const std::vector<std::string>& b) {
        auto tokensA = split(a[0]);
        auto tokensB = split(b[0]);
        if (tokensA.size() < 6 || tokensB.size() < 6) {
            throw std::runtime_error("Alignment first line does not have at least six fields.");
        }
        // Compare by target chromosome (2nd element).
        if (tokensA[1] != tokensB[1])
            return tokensA[1] < tokensB[1];
        // Compare by target start coordinate (3rd element as integer).
        int startA = std::stoi(tokensA[2]);
        int startB = std::stoi(tokensB[2]);
        if (startA != startB)
            return startA < startB;
        // Compare by query chromosome (5th element).
        if (tokensA[4] != tokensB[4])
            return tokensA[4] < tokensB[4];
        // Compare by query start coordinate (6th element as integer).
        int qStartA = std::stoi(tokensA[5]);
        int qStartB = std::stoi(tokensB[5]);
        return qStartA < qStartB;
    });

    // Update the first field of the first line of each alignment block with its new index.
    for (size_t i = 0; i < alignments.size(); i++) {
        auto tokens = split(alignments[i][0]);
        if (tokens.empty()) continue;
        tokens[0] = std::to_string(i);
        alignments[i][0] = join(tokens);
    }

    // Prepare output: print header then sorted alignments.
    std::ostream *outStream = &std::cout;
    std::ofstream outFile;
    if (argc >= 3) {
        outFile.open(argv[2]);
        if (!outFile) {
            std::cerr << "Error opening output file: " << argv[2] << std::endl;
            return 1;
        }
        outStream = &outFile;
    }

    // Write header.
    for (const auto &line : header) {
        (*outStream) << line << "\n";
    }

    // Write each alignment block (each block has exactly 3 lines) with an empty line after each block.
    for (const auto &block : alignments) {
        for (const auto &line : block) {
            (*outStream) << line << "\n";
        }
        (*outStream) << "\n";  // add an empty line between alignments
    }

    if (outFile.is_open()) {
        outFile.close();
    }
    return 0;
}
