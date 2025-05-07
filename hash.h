#ifndef HASH_H
#define HASH_H

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include <algorithm>

typedef std::size_t HASH_INDEX_T;

struct MyStringHash {
    HASH_INDEX_T rValues[5] { 983132572, 1468777056, 552714139, 984953261, 261934300 };
    
    MyStringHash(bool debug = true)
    {
        if(false == debug) {
            generateRValues();
        }
    }

    // hash function entry point (i.e. this is h(k))
    HASH_INDEX_T operator()(const std::string& k) const
    {
        unsigned long long w[5] = {0};  // Store up to 5 base-36 values
        size_t n = k.size();
        
        // Process string in 6-char chunks from the end
        for(size_t i = 0; i < 5; ++i) {
            // Determine substring positions
            size_t endPos = (i * 6 < n) ? (n - i * 6) : 0;
            size_t startPos = (endPos >= 6) ? (endPos - 6) : 0;
            size_t chunkLen = endPos - startPos; // <= 6

            // Build 6-char chunk with 'a' (value 0) padding on the left
            std::string chunk;
            chunk.reserve(6);
            chunk.append(6 - chunkLen, 'a');  // pad with 'a' for zero-value
            chunk.append(k.substr(startPos, chunkLen));

            // Convert chunk from base-36 to integer
            unsigned long long val = 0;
            for(char c : chunk) {
                size_t d = letterDigitToNumber(c);
                val = val * 36 + d;
            }
            w[4 - i] = val;  // last chunk goes into w[4]
        }

        // Combine using rValues
        unsigned long long h = 0;
        for(size_t i = 0; i < 5; ++i) {
            h += w[i] * rValues[i];
        }
        return static_cast<HASH_INDEX_T>(h);
    }

    // convert letter/digit to number 0-35
    HASH_INDEX_T letterDigitToNumber(char letter) const
    {
        letter = std::tolower(static_cast<unsigned char>(letter));
        if(letter >= 'a' && letter <= 'z') return letter - 'a';
        if(letter >= '0' && letter <= '9') return letter - '0' + 26;
        return 0;
    }

    // generate random rValues
    void generateRValues()
    {
        unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());
        std::mt19937 gen(seed);
        for(int i = 0; i < 5; ++i) rValues[i] = gen();
    }
};

#endif // HASH_H
