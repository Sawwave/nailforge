#include "../src/nailforge.hpp"
#include "../argparse/include/argparse/argparse.hpp"
#include <cctype>
#include <string>
#include <cstring>
#include <cstdio>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cmath>

const std::string modelNamePrefix = "TRAIN.";
const bool isTransmark = false;

struct PositionListEntry {
    PositionListEntry(const std::string line) {
        located = false;
        std::istringstream stringStream(line);
        std::string fullModelName;
        stringStream >> fullModelName;

        std::istringstream modelNameStream(fullModelName);
        std::getline(modelNameStream, modelName, '/');

        stringStream >> sequenceName;

        std::string tmp;
        stringStream >> tmp;

        std::string beginString;
        stringStream >> beginString;
        positionFrom = std::stoi(beginString);
        std::string endString;
        stringStream >> endString;
        positionTo = std::stoi(endString);

        if (positionFrom > positionTo) {
            isReverseCompliment = true;
            std::swap(positionFrom, positionTo);
        }
        else {
            isReverseCompliment = false;
        }
    }

    PositionListEntry(const std::string& modelName, const std::string& sequenceName,
        const uint64_t posFrom, const uint64_t posTo) :
        modelName(modelName), sequenceName(sequenceName), positionFrom(posFrom), positionTo(posTo) {
        if (posFrom > posTo) {
            isReverseCompliment = true;
            std::swap(positionFrom, positionTo);
        }
        else {
            isReverseCompliment = false;
        }
    }

    bool operator== (const PositionListEntry& otherEntry)const noexcept {
        const bool hasOverlap = std::max(positionFrom, otherEntry.positionFrom)
            <= std::min(positionTo, otherEntry.positionTo);
        return sequenceName == otherEntry.sequenceName &&
            modelName == otherEntry.modelName &&
            hasOverlap &&
            isReverseCompliment == otherEntry.isReverseCompliment;
    }



    std::string modelName;
    std::string sequenceName;
    uint64_t positionFrom;
    uint64_t positionTo;
    bool isReverseCompliment;
    bool located = false;
};



//Private function prototypes
void parseSearchOptions(argparse::ArgumentParser& parser, std::string& fastaSrc, std::string& awfmiSrc, std::string& hmmSrc,
    std::string& posListSrc, NailForge::SearchParams& params, NailForge::SearchType& searchType, uint8_t& numThreads,
    bool& checkBenchmarkSensitivity, bool& parsableBenchmark, std::string& hitsFileSrc);

void parseCreateOptions(argparse::ArgumentParser& parser, std::string& fastaSrc, std::string& awfmiSrc,
    NailForge::Alphabet& alphabet, uint8_t& suffixArrayCompressionRatio);

bool isPositiveExample(const P7Hmm& phmm, const FastaVector* fastaVector, const uint32_t sequenceIdx);

bool isNegativeExample(const FastaVector* fastaVector, const uint32_t sequenceIdx);

void checkSensitivity(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complimentSeedLiset,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc, const bool parsableBenchmark);

void checkSensitivityTransmark(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complimentSeedLiset, const std::string& posList,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc, const bool parsableBenchmark,
    const uint8_t seedLen, const std::string& hitsFileSrc);


int main(int argc, char** argv) {
    argparse::ArgumentParser parser("nailforgeBenchmark");

    argparse::ArgumentParser createCommand("create");
    createCommand.add_description("create the fm index");
    createCommand.add_argument("-f", "--fasta").required().help("specifies the src for the fast file, for generating the awfm index file");
    createCommand.add_argument("-a", "--awFmIndex").required().help("specifies the src for the awfm index file, for either generating or using the index");
    createCommand.add_argument("-r", "suffixArrayRatio").required().help("suffix array compression ratio to use when generating an awfm index.");
    createCommand.add_argument("-A", "--alphabet").required().help("specifies the alphabet for a fasta. values are d for dna, r for rna, a for amino");

    argparse::ArgumentParser searchCommand("search");
    searchCommand.add_description("searches using a given hmm file and awfmi file");
    searchCommand.add_argument("-f", "--fasta").required().help("specifies the src for the fast file, for generating the awfm index file");
    searchCommand.add_argument("-a", "--awFmIndex").required().help("specifies the src for the awfm index file, for either generating or using the index");
    searchCommand.add_argument("-h", "--hmmFile").required().help("specifies the src for the hmm file, either for generating a model index or for searching");
    searchCommand.add_argument("-l", "--length").required().help("maximum length a diagonal can accumulate to hit the threshold score in order to generate a hit");
    searchCommand.add_argument("-e", "--extension").required().help("length of the extensions used to validate the hits");
    searchCommand.add_argument("-t", "--thresholdScore").required().help("sets the threshold score in bits for the main diagonal hits.");
    searchCommand.add_argument("-x", "--extendedPValue").required().help("P-value for final extended score.");
    searchCommand.add_argument("-n", "--numThreads").help("sets the number of threads when openmp is used.");
    searchCommand.add_argument("-b", "--benchmark").flag().help("benchmark the results for sensitivity");
    searchCommand.add_argument("-B", "--benchmarkData").flag().help("show true positive, false positive results in a easy to parse format");
    searchCommand.add_argument("-p", "--posListSrc").help("file src for the list of real hit positions");
    searchCommand.add_argument("-w", "--extensionWidth").required().help("width of the group of extensions to try");
    searchCommand.add_argument("-H", "--hitsFile").help("file containing the hmmer hits to match against");

    searchCommand.add_argument("-r", "--repeatFilterPartsPerMillion").required().help("float representing maximum number of awfm hits, "
        "in parts per million, before a hit is considered to be repetitive junk");
    searchCommand.add_argument("-T", "--searchType").required().required().help("search type to perform. \n"
        "'s' for standard (amino acid and single strand search),\n "
        "'d' or 'b'. for dual-strand nucleotide search (not to be used with amino acid)\n"
        "'c' for compliment strand only (not to be used with amino acid)");

    parser.add_subparser(createCommand);
    parser.add_subparser(searchCommand);

    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (parser.is_subcommand_used(createCommand)) {
        std::string fastaSrc, awfmSrc;
        uint8_t suffixArrayCompressionRatio;
        NailForge::Alphabet alphabet;

        parseCreateOptions(createCommand, fastaSrc, awfmSrc, alphabet, suffixArrayCompressionRatio);
        std::cout << "Creating fm index..." << std::endl;
        NailForge::ReturnCode nfrc = NailForge::createFmIndex(fastaSrc.c_str(), awfmSrc.c_str(), alphabet, suffixArrayCompressionRatio);

        if (nfrc != NailForge::ReturnCode::Success) {
            std::cerr << "create fm index returned error code " << NailForge::returnCodeDescription(nfrc) << std::endl;
            exit(-1);
        }
        std::cout << "fm index creation finished" << std::endl;
    }

    else if (parser.is_subcommand_used(searchCommand)) {
        std::string fastaSrc, awfmSrc, hmmSrc, posListSrc, hitsFileSrc;
        NailForge::SearchParams params;
        NailForge::SearchType searchType;
        uint8_t numThreads;
        bool checkBenchmarkSensitivity, parsableBenchmark;
        parseSearchOptions(searchCommand, fastaSrc, awfmSrc, hmmSrc, posListSrc, params, searchType,
            numThreads, checkBenchmarkSensitivity, parsableBenchmark, hitsFileSrc);

        std::vector<std::vector<NailForge::AlignmentSeed>> primarySeedList, complimentSeedList;


        // std::cout << "filtering with hmm..." << std::endl;
        auto startTime = std::chrono::high_resolution_clock::now();
        NailForge::ReturnCode nfrc = NailForge::filterWithHmmFile(hmmSrc.c_str(), fastaSrc.c_str(), awfmSrc.c_str(),
            params, searchType, numThreads, primarySeedList, complimentSeedList);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
        const float secondsDuration = (float)duration.count() / 1000000.0f;

        if (checkBenchmarkSensitivity || parsableBenchmark) {
            if (isTransmark) {
                checkSensitivityTransmark(primarySeedList, complimentSeedList, posListSrc, hmmSrc,
                    fastaSrc, parsableBenchmark, params.maximumHitLength, hitsFileSrc);
            }
            else {
                checkSensitivity(primarySeedList, complimentSeedList, hmmSrc, fastaSrc, parsableBenchmark);
            }
        }
        std::cout << "\t" << secondsDuration << std::endl;

    }
    else {
        std::cout << "no command given, please specify createfm or search" << std::endl;
    }

}


void parseCreateOptions(argparse::ArgumentParser& parser, std::string& fastaSrc, std::string& awfmiSrc,
    NailForge::Alphabet& alphabet, uint8_t& suffixArrayCompressionRatio) {
    try {
        fastaSrc = parser.get<std::string>("-f");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get fasta file src from parser" << std::endl;
        std::exit(-1);
    }
    try {
        awfmiSrc = parser.get<std::string>("-a");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get awfmi file src from parser" << std::endl;
        std::exit(-1);
    }
    std::string alphabetStr;
    try {
        alphabetStr = parser.get<std::string>("-A");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get alphabet from parser" << std::endl;
        std::exit(-1);
    }
    if (alphabetStr.length() != 0) {
        switch (std::tolower(alphabetStr[0])) {
        case 'a': alphabet = NailForge::Alphabet::Amino;    break;
        case 'd': alphabet = NailForge::Alphabet::Dna;      break;
        case 'r': alphabet = NailForge::Alphabet::Rna;      break;
        default:
            std::cerr << "Error, alphabet did not start with d, r, or a." << std::endl;
            exit(-1);
        }
    }
    else {
        std::cerr << "alphabet was given an empty string. specify rna, dna, or amino" << std::endl;
        exit(-1);
    }

    try {
        suffixArrayCompressionRatio = std::stoi(parser.get<std::string>("-r"));
    }
    catch (const std::exception& e) {
        std::cerr << "Could not parse suffix array compression ratio" << std::endl;
        exit(-1);
    }
}

void parseSearchOptions(argparse::ArgumentParser& parser, std::string& fastaSrc, std::string& awfmiSrc, std::string& hmmSrc,
    std::string& posListSrc, NailForge::SearchParams& params, NailForge::SearchType& searchType, uint8_t& numThreads,
    bool& checkBenchmarkSensitivity, bool& parsableBenchmark, std::string& hitsFileSrc) {
    try {
        fastaSrc = parser.get<std::string>("-f");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get fasta file src from parser" << std::endl;
        std::exit(-1);
    }
    try {
        awfmiSrc = parser.get<std::string>("-a");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get awfmi file src from parser" << std::endl;
        std::exit(-1);
    }
    try {
        hmmSrc = parser.get<std::string>("-h");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get hmm file src from parser" << std::endl;
        std::exit(-1);
    }
    try {
        hitsFileSrc = parser.get<std::string>("-H");
    }
    catch (const std::exception& e) {
        hitsFileSrc = "";
    }
    try {
        posListSrc = parser.get<std::string>("-p");
    }
    catch (const std::exception& e) {
        //nothing here, not all benchmarks have position list
    }
    try {
        params.maximumHitLength = (uint8_t)std::stoi(parser.get<std::string>("-l"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not get main diag length from parser " << std::endl;
        exit(-1);
    }
    try {
        params.flankExtensionLength = (uint8_t)std::stoi(parser.get<std::string>("-e"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not get extension diag length from parser " << std::endl;
        exit(-1);
    }
    try {
        params.mainDiagonalThresholdScore = std::stof(parser.get<std::string>("-t"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not get main diag threshold from parser " << std::endl;
        exit(-1);
    }
    try {
        params.extensionPValue = std::stof(parser.get<std::string>("-x"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not get extension threshold from parser " << std::endl;
        exit(-1);
    }
    try {
        params.maxSeqHitsPerMillion = std::stof(parser.get<std::string>("-r"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not parse max seq hits per million from parser" << std::endl;
        exit(-1);
    }
    try {
        params.extentionGroupWidth = std::stoi(parser.get<std::string>("-w"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not parse extension group width from parser" << std::endl;
    }
    try {
        numThreads = (uint8_t)std::stoi(parser.get<std::string>("-n"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not get num threads from parser" << std::endl;
        exit(-1);
    }
    try {
        checkBenchmarkSensitivity = parser["-b"] == true;
    }
    catch (const std::exception& e) {
        std::cerr << "could not get benchmark flag from parser" << std::endl;
        exit(-1);
    }
    try {
        parsableBenchmark = parser["-B"] == true;
    }
    catch (const std::exception& e) {
        std::cerr << "could not get parsable benchmark flag from parser" << std::endl;
        exit(-1);
    }
    try {
        std::string searchTypeString = parser.get<std::string>("-T");
        switch (std::tolower(searchTypeString[0])) {
        case 'p': case 'P': searchType = NailForge::SearchType::Standard;           break;
        case 's': case 'S': searchType = NailForge::SearchType::Standard;           break;
        case 'c': case 'C': searchType = NailForge::SearchType::ComplimentStrand;   break;
        case 'b': case 'B': searchType = NailForge::SearchType::DualStrand;         break;
        case 'd': case 'D': searchType = NailForge::SearchType::DualStrand;         break;
        default:
            std::cerr << "search type was not primary, compliment, or dual-strand" << std::endl;
            exit(-1);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "could not parse search type option from parser" << std::endl;
        exit(-1);
    }
}


bool isPositiveExample(const P7Hmm& phmm, const FastaVector* fastaVector, const uint32_t sequenceIdx) {
    size_t headerLen;
    char* headerPtr;
    fastaVectorGetHeader(fastaVector, sequenceIdx, &headerPtr, &headerLen);

    char* modelName = phmm.header.name;
    const auto nameLen = strlen(modelName);

    if (headerLen < nameLen) {
        return false;
    }
    const bool sequenceStartsWithModelName = strncmp(modelName, headerPtr, nameLen) == 0;
    if (sequenceStartsWithModelName) {
        if (headerLen > nameLen) {
            return headerPtr[nameLen] == '/';
        }
        else {
            return true;
        }
    }
    else {
        return false;
    }

}

bool isNegativeExample(const FastaVector* fastaVector, const uint32_t sequenceIdx) {
    size_t headerLen;
    char* headerPtr;
    fastaVectorGetHeader(fastaVector, sequenceIdx, &headerPtr, &headerLen);

    if (headerPtr == NULL) {
        std::cerr << "Error: could not get header from sequence idx " << sequenceIdx << "." << std::endl;
        exit(-1);
    }

    const auto nameLen = strlen("decoy");
    return strncmp("decoy", headerPtr, nameLen) == 0;
}

void checkSensitivity(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complimentSeedLiset,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc, const bool parsableBenchmark) {
    uint64_t numSeeds = 0;
    for (const auto& list : primarySeedList) {
        numSeeds += list.size();
    }

    FastaVector fastaVector;
    FastaVectorReturnCode fvrc = fastaVectorInit(&fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not init fasta vector in sensitivity check" << std::endl;
        exit(-1);
    }
    fvrc = fastaVectorReadFasta(fastaFileSrc.c_str(), &fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not read fasta in sensitivity check" << std::endl;
        exit(-1);
    }
    P7HmmList phmmList;
    P7HmmReturnCode phrc = readP7Hmm(hmmFileSrc.c_str(), &phmmList);
    if (phrc != p7HmmSuccess) {
        std::cerr << "Error: could not read phmm file in sensitivity check" << std::endl;
        exit(-1);
    }


    uint64_t numPositiveExamples = 0;
    uint64_t numNegativeExamples = 0;
    uint64_t numPrimaryTruePositives = 0;
    uint64_t numPrimaryFalsePositives = 0;
    // uint64_t numComplimentTruePositives = 0;
    // uint64_t numComplimentFalsePositives = 0;

    for (uint32_t sequenceIdx = 0; sequenceIdx < fastaVector.metadata.count;sequenceIdx++) {
        if (isNegativeExample(&fastaVector, sequenceIdx)) {
            numNegativeExamples += phmmList.count;;
        }
    }


    for (uint32_t modelIdx = 0;modelIdx < primarySeedList.size(); modelIdx++) {

        for (uint32_t sequenceIdx = 0; sequenceIdx < fastaVector.metadata.count;sequenceIdx++) {
            if (isPositiveExample(phmmList.phmms[modelIdx], &fastaVector, sequenceIdx)) {
                numPositiveExamples++;
                //check to see if this is represented in the seed list
                for (const auto& seed : primarySeedList[modelIdx]) {
                    if (seed.sequenceIdx == sequenceIdx) {
                        numPrimaryTruePositives++;
                        break;
                    }
                }
            }
        }

        //count false positives
        for (const auto& seed : primarySeedList[modelIdx]) {
            if (!isPositiveExample(phmmList.phmms[modelIdx], &fastaVector, seed.sequenceIdx)) {
                numPrimaryFalsePositives++;
            }
        }

    }
    if (numPositiveExamples == 0) {
        std::cout << "error, num pos examples was zero" << std::endl;
    }


    const double primaryTruePositiveRate = static_cast<double>(numPrimaryTruePositives) / static_cast<double>(numPositiveExamples);
    const double primaryFalsePositiveRate = static_cast<double>(numPrimaryFalsePositives) / static_cast<double>(numNegativeExamples);

    if (parsableBenchmark) {
        std::cout << primaryTruePositiveRate  << "\t" << primaryFalsePositiveRate  << "\t"; 
        //<< complimentTruePositiveRate <<   "\t" << complimentFalsePositiveRate;
    }
    else {
        // std::cout << "Primary TP " << numPrimaryTruePositives << " / " << numPositiveExamples <<
        //     " (" << primaryTruePositiveRate * 100 << "%),\tFP " << numPrimaryFalsePositives << " / " << numNegativeExample <<
        //     " (" << primaryFalsePositiveRate * 100 << "%)\n\tCompliment TP " << numComplimentTruePositives << " / " << numPositiveExamples <<
        //     " (" << complimentTruePositiveRate * 100 << "%),\tFP " << numComplimentFalsePositives << " / " << numNegativeExample <<
        //     " (" << complimentFalsePositiveRate * 100 << "%)" << std::endl;
    }



    fastaVectorDealloc(&fastaVector);
    p7HmmListDealloc(&phmmList);
}

void checkSensitivityTransmark(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complimentSeedList, const std::string& posListSrc,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc, const bool parsableBenchmark,
    const uint8_t seedLen, const std::string& hitsFileSrc) {

    //make the position list
    std::ifstream positionStream(posListSrc);
    std::string line;
    std::vector<PositionListEntry> positionList;
    while (std::getline(positionStream, line)) {
        try {
            positionList.emplace_back(line);

        }
        catch (const std::exception& e) {
            std::cout << "ERROR on emplace " << e.what() << std::endl;
        }
    }

    //filter hits out off they're not in the hmmer --max hit list
    std::vector<PositionListEntry> hmmerHits;
    std::ifstream hitsStream(hitsFileSrc);
    while (std::getline(hitsStream, line)) {
        if (line[0] != '#') {
            std::istringstream lineStream(line);
            std::string modelName, sequenceName, modelNameLine;
            uint64_t positionFrom, positionTo;
            std::string tmp;

            lineStream >> sequenceName;
            lineStream >> tmp;
            lineStream >> modelNameLine;
            std::istringstream modelNameStream(modelNameLine);
            std::getline(modelNameStream, tmp, '.');
            modelNameStream >> modelName;

            lineStream >> tmp;
            lineStream >> tmp;
            lineStream >> tmp;
            lineStream >> tmp;
            positionFrom = std::stoi(tmp);
            lineStream >> tmp;
            positionTo = std::stoi(tmp);

            hmmerHits.emplace_back(modelName, sequenceName, positionFrom, positionTo);
        }
    }

    std::vector<PositionListEntry> filteredEntryList;
    for (const auto& position : positionList) {
        for (const auto& hmmerHit : hmmerHits) {
            if (position == hmmerHit) {
                filteredEntryList.push_back(position);
                break;
            }
        }
    }


    FastaVector fastaVector;
    FastaVectorReturnCode fvrc = fastaVectorInit(&fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not init fasta vector in sensitivity check" << std::endl;
        exit(-1);
    }
    fvrc = fastaVectorReadFasta(fastaFileSrc.c_str(), &fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not read fasta in sensitivity check" << std::endl;
        exit(-1);
    }
    P7HmmList phmmList;
    P7HmmReturnCode phrc = readP7Hmm(hmmFileSrc.c_str(), &phmmList);
    if (phrc != p7HmmSuccess) {
        std::cerr << "Error: could not read phmm file in sensitivity check" << std::endl;
        exit(-1);
    }


    uint64_t numPrimaryTruePositives = 0;
    uint64_t numComplimentTruePositives = 0;
    uint64_t numPrimaryPositiveExamples = 0;
    uint64_t numComplimentPositiveExamples = 0;
    uint64_t numPrimarySeedsFromAllModels = 0;
    uint64_t numComplimentSeedsFromAllModels = 0;
    uint64_t numImpliedSequences = 0;

    for (uint64_t modelIdx = 0; modelIdx < phmmList.count;modelIdx++) {
        numImpliedSequences += (fastaVector.sequence.count / phmmList.phmms[modelIdx].header.modelLength);
        numPrimarySeedsFromAllModels += primarySeedList[modelIdx].size();
        numComplimentSeedsFromAllModels += complimentSeedList[modelIdx].size();

        const uint32_t modelLength = phmmList.phmms[modelIdx].header.modelLength;
        const char* modelName = phmmList.phmms[modelIdx].header.name;
        std::vector<PositionListEntry> thisModelsPrimaryEntries;
        std::vector<PositionListEntry> thisModelComplimentEntries;

        //get only the real hits that correspond to this model.
        std::copy_if(filteredEntryList.begin(), filteredEntryList.end(), std::back_inserter(thisModelsPrimaryEntries),
            [&](PositionListEntry& p) {return std::strcmp((modelNamePrefix + p.modelName).c_str(), modelName) == 0 && !p.isReverseCompliment;});
        std::copy_if(filteredEntryList.begin(), filteredEntryList.end(), std::back_inserter(thisModelComplimentEntries),
            [&](PositionListEntry& p) {return std::strcmp((modelNamePrefix + p.modelName).c_str(), modelName) == 0 && p.isReverseCompliment;});

        numPrimaryPositiveExamples += thisModelsPrimaryEntries.size();
        numComplimentPositiveExamples += thisModelComplimentEntries.size();

        for (auto& positionEntry : thisModelsPrimaryEntries) {

            for (const auto seed : primarySeedList[modelIdx]) {
                char* seqHeader;
                uint64_t headerLen;
                fastaVectorGetHeader(&fastaVector, seed.sequenceIdx, &seqHeader, &headerLen);
                if (strcmp(positionEntry.sequenceName.c_str(), seqHeader) == 0) {
                    //check the positions
                    const uint64_t seedEndPos = seed.sequencePosition + seedLen;
                    bool hasOverlap = std::max(seed.sequencePosition, positionEntry.positionFrom) <=
                        std::min(seedEndPos, positionEntry.positionTo);
                    if (hasOverlap) {
                        positionEntry.located = true;
                    }
                }
            }
        }

        for (auto& positionEntry : thisModelComplimentEntries) {

            for (const auto seed : complimentSeedList[modelIdx]) {
                char* seqHeader;
                uint64_t headerLen;
                fastaVectorGetHeader(&fastaVector, seed.sequenceIdx, &seqHeader, &headerLen);
                if (strcmp(positionEntry.sequenceName.c_str(), seqHeader) == 0) {
                    //check the positions
                    const uint64_t seedEndPos = seed.sequencePosition + seedLen;
                    bool hasOverlap = std::max(seed.sequencePosition, positionEntry.positionFrom) <=
                        std::min(seedEndPos, positionEntry.positionTo);
                    if (hasOverlap) {
                        positionEntry.located = true;
                    }
                }
            }
        }
        const auto thisModelPrimaryTruePositives = std::count_if(thisModelsPrimaryEntries.begin(), thisModelsPrimaryEntries.end(),
            [](PositionListEntry& p) {return p.located;});
        numPrimaryTruePositives += thisModelPrimaryTruePositives;
        const auto thisModelComplimentTruePositives = std::count_if(thisModelComplimentEntries.begin(), thisModelComplimentEntries.end(),
            [](PositionListEntry& p) {return p.located;});
        numComplimentTruePositives += thisModelComplimentTruePositives;
    }

    const uint64_t numPrimaryFalsePositives = numPrimarySeedsFromAllModels - numPrimaryTruePositives;
    const uint64_t numComplimentFalsePositives = numComplimentSeedsFromAllModels - numComplimentTruePositives;
    const double primaryTpRate = static_cast<double>(numPrimaryTruePositives) / static_cast<double>(numPrimaryPositiveExamples);
    const double complimentTpRate = static_cast<double>(numComplimentTruePositives) / static_cast<double>(numComplimentPositiveExamples);
    const double primaryFpRate = static_cast<double>(numPrimaryFalsePositives) / static_cast<double>(numImpliedSequences - numPrimaryTruePositives);
    const double complimentFpRate = static_cast<double>(numComplimentFalsePositives) / static_cast<double>(numImpliedSequences - numComplimentTruePositives);

    if (parsableBenchmark) {
        std::cout << primaryTpRate << "\t" << primaryFpRate << "\t" << complimentTpRate <<
            "\t" << complimentFpRate;
    }
    else {
        std::cout << "Primary TP rate" << primaryTpRate * 100.0f << "%\t Primary FP rate " <<
            primaryFpRate * 100.0f << "%\tCompliment TP rate " << complimentTpRate * 100.0f <<
            "%\t Compliment FP rate" << complimentFpRate * 100.0f << "%" << std::endl;
    }

    fastaVectorDealloc(&fastaVector);
    p7HmmListDealloc(&phmmList);

}