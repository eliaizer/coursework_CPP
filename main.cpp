#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <array>
#include <optional>
#include <algorithm>
#include <numeric>
#include <map>
#include <chrono>
#include <filesystem>
#include <string>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

namespace util {
    template<class T>
    bool read(std::istream &in, T &value) {
        return static_cast<bool>(in.read(reinterpret_cast<char*>(&value), sizeof(T)));
    }

    template<class T>
    void write(std::ostream &out, const T &value) {
        out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }
}

class HuffmanCodec {
    struct Node {
        uint32_t freq;
        int16_t symbol;
        Node *left{}, *right{};
        Node(uint32_t f, int16_t s, Node *l=nullptr, Node *r=nullptr) : freq(f), symbol(s), left(l), right(r) {}
    };
    struct cmp {
        bool operator()(const Node* a, const Node* b) const { return a->freq > b->freq; }
    };
    using CodeTable = std::array<std::pair<uint64_t, uint8_t>, 256>;

public:
    void encodeFile(const std::string &inPath, const std::string &outPath);
    void decodeFile(const std::string &inPath, const std::string &outPath);

private:
    std::array<uint32_t,256> buildFrequency(const std::vector<uint8_t> &data);
    Node* buildTree(const std::array<uint32_t,256> &freq);
    void buildCodeTable(Node *n, uint64_t code, uint8_t len, CodeTable &table);
    void deleteTree(Node *n);
};

void HuffmanCodec::encodeFile(const std::string &inPath, const std::string &outPath) {
    std::ifstream fin(inPath, std::ios::binary);
    std::vector<uint8_t> data((std::istreambuf_iterator<char>(fin)), {});
    auto freq = buildFrequency(data);
    Node *root = buildTree(freq);
    CodeTable table{};
    buildCodeTable(root, 0, 0, table);
    std::ofstream fout(outPath, std::ios::binary);
    for (uint32_t f : freq) util::write(fout, f);
    uint64_t bitBuf = 0; int bitsStored = 0;
    auto flushBits = [&]() {
        while (bitsStored >= 8) {
            fout.put(static_cast<char>(bitBuf & 0xFF));
            bitBuf >>= 8; bitsStored -= 8;
        }
    };
    for (uint8_t b : data) {
        auto [code, len] = table[b];
        bitBuf |= (code << bitsStored);
        bitsStored += len;
        flushBits();
    }
    uint8_t remaining = bitsStored;
    flushBits();
    if (bitsStored) fout.put(static_cast<char>(bitBuf & 0xFF));
    fout.put(static_cast<char>(remaining));
    deleteTree(root);
}

void HuffmanCodec::decodeFile(const std::string &inPath, const std::string &outPath) {
    std::ifstream fin(inPath, std::ios::binary);
    std::array<uint32_t,256> freq{};
    for (uint32_t &f : freq) util::read(fin, f);
    Node *root = buildTree(freq);
    std::vector<uint8_t> body((std::istreambuf_iterator<char>(fin)), {});
    uint8_t remainingBits = body.back(); body.pop_back();
    std::ofstream fout(outPath, std::ios::binary);
    Node *cur = root;
    size_t totalBytes = std::accumulate(freq.begin(), freq.end(), size_t{0});
    size_t written = 0;
    for (size_t i = 0; i < body.size(); ++i) {
        uint8_t byte = body[i];
        int bitCount = (i + 1 == body.size()) ? remainingBits : 8;
        for (int bit = 0; bit < bitCount; ++bit) {
            cur = (byte & 1) ? cur->right : cur->left;
            byte >>= 1;
            if (!cur->left && !cur->right) {
                fout.put(static_cast<char>(cur->symbol));
                if (++written == totalBytes) break;
                cur = root;
            }
        }
    }
    deleteTree(root);
}

std::array<uint32_t,256> HuffmanCodec::buildFrequency(const std::vector<uint8_t> &data) {
    std::array<uint32_t,256> freq{};
    for (uint8_t b : data) ++freq[b];
    return freq;
}

HuffmanCodec::Node* HuffmanCodec::buildTree(const std::array<uint32_t,256> &freq) {
    std::priority_queue<Node*, std::vector<Node*>, cmp> pq;
    for (int i = 0; i < 256; ++i) if (freq[i]) pq.push(new Node(freq[i], i));
    if (pq.empty()) pq.push(new Node(1, 0));
    while (pq.size() > 1) {
        Node *a = pq.top(); pq.pop();
        Node *b = pq.top(); pq.pop();
        pq.push(new Node(a->freq + b->freq, -1, a, b));
    }
    return pq.top();
}

void HuffmanCodec::buildCodeTable(Node *n, uint64_t code, uint8_t len, CodeTable &table) {
    if (!n->left && !n->right) { table[n->symbol] = {code, std::max<uint8_t>(1, len)}; return; }
    buildCodeTable(n->left,  code, len + 1, table);
    buildCodeTable(n->right, code | (1ULL << len), len + 1, table);
}

void HuffmanCodec::deleteTree(Node *n) {
    if (!n) return;
    deleteTree(n->left);
    deleteTree(n->right);
    delete n;
}

class LZ77Codec {
    static constexpr size_t WINDOW_SIZE = 4096, LOOKAHEAD = 18;
    struct Token { uint16_t offset; uint8_t length; uint8_t next; };
public:
    void encodeFile(const std::string &inPath, const std::string &outPath);
    void decodeFile(const std::string &inPath, const std::string &outPath);
private:
    Token findLongestMatch(const std::vector<uint8_t> &data, size_t pos);
};

void LZ77Codec::encodeFile(const std::string &inPath, const std::string &outPath) {
    std::ifstream fin(inPath, std::ios::binary);
    std::vector<uint8_t> data((std::istreambuf_iterator<char>(fin)), {});
    std::ofstream fout(outPath, std::ios::binary);
    size_t pos = 0;
    while (pos < data.size()) {
        auto token = findLongestMatch(data, pos);
        util::write(fout, token.offset);
        fout.put(static_cast<char>(token.length));
        fout.put(static_cast<char>(token.next));
        pos += token.length + 1;
    }
    util::write<uint16_t>(fout, 0); fout.put(0); fout.put(0);
}

void LZ77Codec::decodeFile(const std::string &inPath, const std::string &outPath) {
    std::ifstream fin(inPath, std::ios::binary);
    std::vector<uint8_t> outBuf;
    while (true) {
        uint16_t offset; if (!util::read(fin, offset)) break;
        uint8_t length;  fin.read(reinterpret_cast<char*>(&length), 1);
        uint8_t next;    fin.read(reinterpret_cast<char*>(&next), 1);
        if (offset==0 && length==0) break;
        if (offset==0 || length==0) { outBuf.push_back(next); continue; }
        size_t start = outBuf.size() - offset;
        for (size_t i = 0; i < length; ++i)
            outBuf.push_back(outBuf[start + i]);
        outBuf.push_back(next);
    }
    std::ofstream fout(outPath, std::ios::binary);
    fout.write(reinterpret_cast<const char*>(outBuf.data()), outBuf.size());
}

LZ77Codec::Token LZ77Codec::findLongestMatch(const std::vector<uint8_t> &data, size_t pos) {
    size_t bestLen = 0, bestOff = 0;
    size_t windowStart = (pos > WINDOW_SIZE) ? pos - WINDOW_SIZE : 0;
    for (size_t i = windowStart; i < pos; ++i) {
        size_t len = 0;
        while (len < LOOKAHEAD && pos + len < data.size() && data[i + len] == data[pos + len]) ++len;
        if (len > bestLen && len >= 3) { bestLen = len; bestOff = pos - i; if (len == LOOKAHEAD) break; }
    }
    if (bestLen >= 3) return {static_cast<uint16_t>(bestOff), static_cast<uint8_t>(bestLen), data[pos + bestLen]};
    return {0, 0, data[pos]};
}

class LZWCodec {
public:
    void encodeFile(const std::string &inPath, const std::string &outPath);
    void decodeFile(const std::string &inPath, const std::string &outPath);
};

void LZWCodec::encodeFile(const std::string &inPath, const std::string &outPath) {
    std::ifstream fin(inPath, std::ios::binary);
    std::vector<uint8_t> data((std::istreambuf_iterator<char>(fin)), {});
    std::ofstream fout(outPath, std::ios::binary);
    std::unordered_map<std::string, int> dict;
    for (int i = 0; i < 256; ++i) dict[std::string(1, i)] = i;
    std::string w;
    int code = 256;
    for (uint8_t c : data) {
        std::string wc = w + static_cast<char>(c);
        if (dict.count(wc)) w = wc;
        else {
            util::write(fout, static_cast<uint16_t>(dict[w]));
            dict[wc] = code++;
            w = std::string(1, c);
        }
    }
    if (!w.empty()) util::write(fout, static_cast<uint16_t>(dict[w]));
}

void LZWCodec::decodeFile(const std::string &inPath, const std::string &outPath) {
    std::ifstream fin(inPath, std::ios::binary);
    std::ofstream fout(outPath, std::ios::binary);
    std::map<int, std::string> dict;
    for (int i = 0; i < 256; ++i) dict[i] = std::string(1, i);
    uint16_t prevCode, currCode;
    util::read(fin, prevCode);
    fout << dict[prevCode];
    std::string w = dict[prevCode];
    int code = 256;
    while (util::read(fin, currCode)) {
        std::string entry;
        if (dict.count(currCode)) entry = dict[currCode];
        else entry = w + w[0];
        fout << entry;
        dict[code++] = w + entry[0];
        w = entry;
    }
}

std::string buildOutputFileName(const std::string& inputFile, const std::string& algoSuffix) {
    std::filesystem::path p(inputFile);
    auto dir = p.parent_path();
    auto stem = p.stem();
    auto ext = p.extension();
    std::string outFileName = stem.string() + "_" + algoSuffix + ext.string();
    return (dir / outFileName).string();
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 5) {
        std::cerr << "Проверьте аргументы: ./compress <huff|lz77|lzw> <encode|decode> <input>\n";
        return 1;
    }
    std::string algo = argv[1];
    std::string mode = argv[2];
    std::string inFile = argv[3];
    std::string outFile;

    std::string suffix;
    if (algo == "huff") suffix = "huff";
    else if (algo == "lz77") suffix = "lz77";
    else if (algo == "lzw") suffix = "lzw";
    else return 1;
    outFile = buildOutputFileName(inFile, suffix);

    try {
        if (algo == "huff") {
            HuffmanCodec c;
            if (mode == "encode")
                c.encodeFile(inFile, outFile);
            else
                c.decodeFile(inFile, outFile);
        } else if (algo == "lz77") {
            LZ77Codec c;
            if (mode == "encode")
                c.encodeFile(inFile, outFile);
            else
                c.decodeFile(inFile, outFile);
        } else if (algo == "lzw") {
            LZWCodec c;
            if (mode == "encode")
                c.encodeFile(inFile, outFile);
            else
                c.decodeFile(inFile, outFile);
        } else {
            return 1;
        }
    } catch (...) {
        return 2;
    }
    return 0;
}
