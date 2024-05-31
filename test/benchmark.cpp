#include "../src/nailforge.hpp"
#include "../argparse/include/argparse/argparse.hpp"
#include "../src/PhmmProcessor/PhmmProcessor.hpp"
#include "../src/Alphabet/LetterConversion.hpp"
#include <cctype>
#include <string>
#include <cstring>
#include <cstdio>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cmath>
#include <optional>

using std::vector;
using std::string;
using std::unique_ptr;

const std::string modelNamePrefix = "TRAIN.";

struct PositionListEntry {
    PositionListEntry(const std::string line) {
        located = false;
        std::istringstream stringStream(line);
        std::string fullModelName;
        std::string tmp;
        std::string beginString;
        std::string endString;

        stringStream >> fullModelName;

        std::istringstream modelNameStream(fullModelName);
        std::getline(modelNameStream, modelName, '/');

        stringStream >> sequenceName;

        stringStream >> tmp;

        stringStream >> beginString;
        positionFrom = std::stoi(beginString);
        stringStream >> endString;
        positionTo = std::stoi(endString);

        if (positionFrom > positionTo) {
            isReverseComplement = true;
            std::swap(positionFrom, positionTo);
        }
        else {
            isReverseComplement = false;
        }
    }

    PositionListEntry(const std::string& modelName, const std::string& sequenceName,
        const uint64_t posFrom, const uint64_t posTo) :
        modelName(modelName), sequenceName(sequenceName), positionFrom(posFrom), positionTo(posTo) {
        if (posFrom > posTo) {
            isReverseComplement = true;
            std::swap(positionFrom, positionTo);
        }
        else {
            isReverseComplement = false;
        }
    }

    bool operator== (const PositionListEntry& otherEntry)const noexcept {
        const bool hasOverlap = std::max(positionFrom, otherEntry.positionFrom)
            <= std::min(positionTo, otherEntry.positionTo);
        return sequenceName == otherEntry.sequenceName &&
            modelName == otherEntry.modelName &&
            hasOverlap &&
            isReverseComplement == otherEntry.isReverseComplement;
    }

    std::string modelName;
    std::string sequenceName;
    uint64_t positionFrom;
    uint64_t positionTo;
    bool isReverseComplement;
    bool located = false;
};



//Private function prototypes
void parseSearchOptions(argparse::ArgumentParser& parser, std::string& fastaSrc, std::string& awfmiSrc, std::string& hmmSrc,
    std::string& posListSrc, NailForge::SearchParams& params, NailForge::SearchType& searchType, uint8_t& numThreads,
    NailForge::Alphabet& searchAlphabet);

void parseCreateOptions(argparse::ArgumentParser& parser, std::string& fastaSrc, std::string& awfmiSrc,
    NailForge::Alphabet& alphabet, uint8_t& suffixArrayCompressionRatio);

bool isPositiveExample(const P7Hmm& phmm, const FastaVector* fastaVector, const uint32_t sequenceIdx);

bool isAminoNegativeExample(const FastaVector* fastaVector, const uint32_t sequenceIdx);

void checkAminoSensitivity(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc, const std::string& posListSrc);

void checkDnaSensitivity(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complementSeedLiset, const std::string& posList,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc);

unique_ptr<vector<vector<PositionListEntry>>> getDnaPositionList(const std::string& positionListFileSrc, const P7HmmList& phmmList, const FastaVector& fastaVector);

std::vector<std::vector<PositionListEntry>> filterDnaPositionLists(const std::vector<std::vector<PositionListEntry>> hmmerHits,
    const std::vector<std::vector<PositionListEntry>> officialPositionList);

void writeResultsToFile(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complementSeedList, const std::string& hmmSrc, const std::string& fastaSrc);

void writeHitsToFile(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complementSeedList, const std::string& fastaSrc,
    const std::string hmmSrc, const std::string& outputFileSrc);

void writeAminoScoreList(std::vector<std::vector<NailForge::AlignmentSeed>>& seedList, const std::string& fastaSrc,
    const std::string hmmSrc);

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
    searchCommand.add_argument("-p", "--posListSrc").help("file src for the list of real hit positions");
    searchCommand.add_argument("-A", "--alphabet").required().help("sets the alphabet, necessary for determining how to parse the TP/FP rates");
    searchCommand.add_argument("-r", "--repeatFilterPartsPerMillion").required().help("float representing maximum number of awfm hits, "
        "in parts per million, before a hit is considered to be repetitive junk");
    searchCommand.add_argument("-T", "--searchType").required().help("search type to perform. \n"
        "'s' for standard (amino acid and single strand search),\n "
        "'d' or 'b'. for dual-strand nucleotide search (not to be used with amino acid)\n"
        "'c' for complement strand only (not to be used with amino acid)");

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
        std::string fastaSrc, awfmSrc, hmmSrc, posListSrc;
        NailForge::SearchParams params;
        NailForge::SearchType searchType;
        uint8_t numThreads;
        NailForge::Alphabet searchAlphabet;

        parseSearchOptions(searchCommand, fastaSrc, awfmSrc, hmmSrc, posListSrc, params, searchType,
            numThreads, searchAlphabet);

        std::vector<std::vector<NailForge::AlignmentSeed>> primarySeedList, complementSeedList;

        auto startTime = std::chrono::high_resolution_clock::now();
        float searchTimeDuration = 0.0f;
        NailForge::ReturnCode nfrc = NailForge::filterWithHmmFile(hmmSrc.c_str(), fastaSrc.c_str(), awfmSrc.c_str(),
            params, searchType, numThreads, primarySeedList, complementSeedList, searchTimeDuration);
        if (nfrc != NailForge::ReturnCode::Success) {
            std::cerr << "unexpected error from filterWithHmm: " << NailForge::returnCodeDescription(nfrc) << std::endl;
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto fullRuntimeDuration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
        const float fullRuntimeDurationSeconds = (float)fullRuntimeDuration.count() / 1000000.0f;


        std::cout << "done writing" << std::endl;
        if (searchAlphabet == NailForge::Alphabet::Amino) {
            checkAminoSensitivity(primarySeedList, hmmSrc, fastaSrc, posListSrc);
        }
        else {
            checkDnaSensitivity(primarySeedList, complementSeedList, posListSrc, hmmSrc, fastaSrc);
        }


        std::cout << "\ttime (search): \t" << searchTimeDuration << "\t time(full)\t" << fullRuntimeDurationSeconds << std::endl;
        
        std::cout << "writing hits to file" << std::endl;
        writeHitsToFile(primarySeedList, complementSeedList, fastaSrc, hmmSrc, "benchmarkHits.txt");
        writeAminoScoreList(primarySeedList, fastaSrc, hmmSrc);
        writeResultsToFile(primarySeedList, complementSeedList, hmmSrc, fastaSrc);

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
    NailForge::Alphabet& searchAlphabet) {
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
        posListSrc = parser.get<std::string>("-p");
    }
    catch (const std::exception& e) {
        //nothing here, not all benchmarks have position list
    }
    try {
        params.maximumHitLength = (uint32_t)std::stoi(parser.get<std::string>("-l"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not get main diag length from parser " << std::endl;
        exit(-1);
    }
    try {
        params.flankExtensionLength = (uint32_t)std::stoi(parser.get<std::string>("-e"));
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
        numThreads = (uint8_t)std::stoi(parser.get<std::string>("-n"));
    }
    catch (const std::exception& e) {
        std::cerr << "could not get num threads from parser" << std::endl;
        exit(-1);
    }
    try {
        std::string searchTypeString = parser.get<std::string>("-T");
        switch (std::tolower(searchTypeString[0])) {
        case 'p': case 'P': searchType = NailForge::SearchType::Standard;           break;
        case 's': case 'S': searchType = NailForge::SearchType::Standard;           break;
        case 'c': case 'C': searchType = NailForge::SearchType::ComplementStrand;   break;
        case 'b': case 'B': searchType = NailForge::SearchType::DualStrand;         break;
        case 'd': case 'D': searchType = NailForge::SearchType::DualStrand;         break;
        default:
            std::cerr << "search type was not primary, complement, or dual-strand" << std::endl;
            exit(-1);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "could not parse search type option from parser" << std::endl;
        exit(-1);
    }
    try {
        std::string searchAlphabetString = parser.get<std::string>("-A");
        switch (std::tolower(searchAlphabetString[0])) {
        case 'd':  searchAlphabet = NailForge::Alphabet::Dna;  break;
        case 'r': searchAlphabet = NailForge::Alphabet::Rna;   break;
        case 'a': searchAlphabet = NailForge::Alphabet::Amino; break;
        default:
            std::cerr << "search alphabet must be dna, rna, or amino." << std::endl;
            exit(-1);
        }
        if (searchAlphabet == NailForge::Alphabet::Amino && searchType != NailForge::SearchType::Standard) {
            std::cerr << "error, only Standard search type is supported for amino acid search" << std::endl;
            exit(-1);
        }
    }
    catch (const std::exception& e) {
        std::cerr << "could not parse alphabet option" << std::endl;
        exit(-1);
    }
}


bool isPositiveExample(const P7Hmm& phmm, const FastaVector* fastaVector, const uint32_t sequenceIdx) {
    size_t headerLen;
    char* headerPtr;
    fastaVectorGetHeader(fastaVector, sequenceIdx, &headerPtr, &headerLen);

    char* modelName = phmm.header.name;

    bool foundModelNameInHeader = std::string(modelName) == std::string(headerPtr).substr(0, std::string(headerPtr).find("/"));
    return foundModelNameInHeader;
    // const auto nameLen = strlen(modelName);

    // if (headerLen < nameLen) {
    //     return false;
    // }
    // const bool sequenceStartsWithModelName = strncmp(modelName, headerPtr, nameLen) == 0;
    // if (sequenceStartsWithModelName) {
    //     if (headerLen > nameLen) {
    //         return headerPtr[nameLen] == '/';
    //     }
    //     else {
    //         return true;
    //     }
    // }
    // else {
    //     return false;
    // }

}


bool isAminoNegativeExample(const FastaVector* fastaVector, const uint32_t sequenceIdx) {
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

void checkAminoSensitivity(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc, const std::string& posListSrc) {
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
    uint64_t numTruePositives = 0;
    uint64_t numFalsePositives = 0;

    //find the number of negative examples
    for (uint32_t sequenceIdx = 0; sequenceIdx < fastaVector.metadata.count;sequenceIdx++) {
        if (isAminoNegativeExample(&fastaVector, sequenceIdx)) {
            numNegativeExamples += phmmList.count;
        }
    }

    //find the false positives
    for (uint32_t modelIdx = 0; modelIdx < primarySeedList.size(); modelIdx++) {
        //count false positives
        for (const auto& seed : primarySeedList[modelIdx]) {
            if (isAminoNegativeExample(&fastaVector, seed.sequenceIdx)) {
                numFalsePositives++;
            }
        }
    }

    //find the positive examples and true positives
    std::string line;

    std::ifstream positionListFile(posListSrc);
    while (std::getline(positionListFile, line)) {
        if (line[0] != '#') {
            numPositiveExamples++;
            size_t stringPosition = line.find('/');
            stringPosition++;
            stringPosition = line.find('/', stringPosition);
            const std::string sequenceName = line.substr(0, stringPosition + 1);

            for (uint32_t modelIdx = 0; modelIdx < primarySeedList.size(); modelIdx++) {
                char* modelName = phmmList.phmms[modelIdx].header.name;

                if (sequenceName.find(std::string(modelName)) != 0) {
                    continue;
                }

                //find the index of the matching sequence
                for (uint32_t sequenceIdx = 0; sequenceIdx < fastaVector.metadata.count;sequenceIdx++) {
                    char* seqHeader;
                    size_t seqHeaderLen;
                    fastaVectorGetHeader(&fastaVector, sequenceIdx, &seqHeader, &seqHeaderLen);
                    if (std::string(seqHeader).find(sequenceName) != std::string::npos) {
                        for (const NailForge::AlignmentSeed& primarySeed : primarySeedList[modelIdx]) {
                            if (primarySeed.sequenceIdx == sequenceIdx) {
                                numTruePositives++;
                                break;
                            }
                        }
                        break;
                    }

                }
                break;//we've identified the model, we can skip the rest
            }
        }
    }
    if (numPositiveExamples == 0) {
        std::cout << "error, num pos examples was zero" << std::endl;
    }

    size_t numPrimarySeeds = 0;
    for (const auto& primaryHit : primarySeedList) {
        numPrimarySeeds += primaryHit.size();
    }
    const double primaryTruePositiveRate = static_cast<double>(numTruePositives) / static_cast<double>(numPositiveExamples);
    const double primaryFalsePositiveRate = static_cast<double>(numFalsePositives) / static_cast<double>(numNegativeExamples);

    std::cout << "#primary hits: " << numPrimarySeeds << "\t";
    std::cout << "TP:\t" << numTruePositives << "/" << numPositiveExamples << "\t(" << primaryTruePositiveRate << ")\tFP:\t" <<
        numFalsePositives << "/" << numNegativeExamples << "\t(" << primaryFalsePositiveRate << ")";
    fastaVectorDealloc(&fastaVector);
    p7HmmListDealloc(&phmmList);
}

void checkDnaSensitivity(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complementSeedList, const std::string& posListSrc,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc) {

    //make the position list
    std::ifstream positionStream(posListSrc);
    std::string line;


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


    const unique_ptr<vector<vector<PositionListEntry>>> positionList = getDnaPositionList(posListSrc, phmmList, fastaVector);


    double numImpliedSequences = 0.0f;
    uint64_t numPrimaryTruePositives = 0;
    uint64_t numComplementTruePositives = 0;
    uint64_t numPrimaryFalsePositives = 0;
    uint64_t numComplementFalsePositives = 0;
    uint64_t numPrimaryPositiveExamples = 0;
    uint64_t numComplementPositiveExamples = 0;
    uint64_t numNegativeExamples = 0;


    uint64_t numPosNotFound = 0;
    //primary strand count true positives
    for (uint64_t modelIdx = 0; modelIdx < phmmList.count;modelIdx++) {
        //this seems hackish and difficult to explain, can we do better?
        numImpliedSequences += static_cast<double>(fastaVector.sequence.count) /
            static_cast<double>(phmmList.phmms[modelIdx].header.modelLength);

        const auto& posList = *positionList;
        const auto& modPosList = posList[modelIdx];
        for (const auto& position : modPosList) {
            bool foundPosition = false;
            if (position.isReverseComplement) {
                continue;
            }
            for (const auto& nailHit : primarySeedList[modelIdx]) {
                char* seqName;
                uint64_t seqNameLen;
                fastaVectorGetHeader(&fastaVector, nailHit.sequenceIdx, &seqName, &seqNameLen);
                bool sequencesMatch = std::strcmp(seqName, position.sequenceName.c_str()) == 0;
                bool hasOverlap = std::max((int64_t)position.positionFrom, std::max((int64_t)nailHit.sequencePosition-128,(int64_t)0 )) <=
                    std::min((int64_t)position.positionTo, (int64_t)nailHit.sequencePosition + 128);

                // nailHit.sequencePosition >= position.positionFrom && nailHit.sequencePosition <= position.positionTo;
                bool seedMatchesActualPosition = sequencesMatch && hasOverlap;
                if (seedMatchesActualPosition) {
                    numPrimaryTruePositives++;
                    foundPosition = true;
                    break;
                }
            }
            if(!foundPosition){
                numPosNotFound++;
                std::cout << "couldn't find position :"<< position.sequenceName<<"\t"<<position.modelName<<"\t"<< position.positionFrom << "\t"<<position.positionTo<< "\t"<<position.isReverseComplement << std::endl; 
            }

        }
    }
    std::cout << "couldn't find "<< numPosNotFound<<" positions"<<std::endl;

    //complementStrand count true positives
    for (uint32_t modelIdx = 0; modelIdx < phmmList.count;modelIdx++) {
        auto& posListForModel = (*positionList)[modelIdx];
        for (auto& position : posListForModel) {
            if (!position.isReverseComplement) {
                continue;
            }
            for (const auto& nailHit : complementSeedList[modelIdx]) {
                char* seqName;
                uint64_t seqNameLen;
                fastaVectorGetHeader(&fastaVector, nailHit.sequenceIdx, &seqName, &seqNameLen);
                bool sequencesMatch = std::strcmp(seqName, position.sequenceName.c_str()) == 0;
                bool hasOverlap = nailHit.sequencePosition >= position.positionFrom && nailHit.sequencePosition <= position.positionTo;
                bool seedMatchesActualPosition = sequencesMatch && hasOverlap;
                if (seedMatchesActualPosition) {
                    numComplementTruePositives++;
                    break;
                }
            }
        }
    }


    //count primary false positives
    for (uint32_t modelIdx = 0; modelIdx < phmmList.count;modelIdx++) {
        for (const auto& nailHit : primarySeedList[modelIdx]) {
            bool foundCorrespondingHit = false;
            bool matchesAgainstAnotherModel = false;
            for (const auto& position : (*positionList)[modelIdx]) {
                char* seqName;
                uint64_t seqNameLen;
                fastaVectorGetHeader(&fastaVector, nailHit.sequenceIdx, &seqName, &seqNameLen);
                bool sequencesMatch = std::strcmp(seqName, position.sequenceName.c_str()) == 0;
                bool hasOverlap = nailHit.sequencePosition >= position.positionFrom && nailHit.sequencePosition <= position.positionTo;
                if (sequencesMatch && hasOverlap && !position.isReverseComplement) {
                    foundCorrespondingHit = true;
                    break;
                }
                else if (sequencesMatch && hasOverlap) {
                    matchesAgainstAnotherModel = true;
                }
            }
            if (!foundCorrespondingHit && !matchesAgainstAnotherModel) {
                numPrimaryFalsePositives++;
            }
        }
    }

    //count complement false positives
    for (uint32_t modelIdx = 0; modelIdx < phmmList.count;modelIdx++) {
        for (const auto& nailHit : complementSeedList[modelIdx]) {
            bool foundCorrespondingHit = false;
            for (const auto& position : (*positionList)[modelIdx]) {
                char* seqName;
                uint64_t seqNameLen;
                fastaVectorGetHeader(&fastaVector, nailHit.sequenceIdx, &seqName, &seqNameLen);
                bool sequencesMatch = std::strcmp(seqName, position.sequenceName.c_str()) == 0;
                bool hasOverlap = nailHit.sequencePosition >= position.positionFrom && nailHit.sequencePosition <= position.positionTo;
                if (sequencesMatch && hasOverlap && position.isReverseComplement) {
                    foundCorrespondingHit = true;
                    break;
                }
            }
            if (!foundCorrespondingHit) {
                numComplementFalsePositives++;
            }
        }
    }

    //count primary/complement positive examples
    for (uint32_t modelIdx = 0; modelIdx < phmmList.count;modelIdx++) {
        for (const auto& position : (*positionList)[modelIdx]) {
            if (position.isReverseComplement) {
                numComplementPositiveExamples++;
            }
            else {
                numPrimaryPositiveExamples++;
            }
        }
    }


    //count number of negative examples
    for (uint64_t modelIdx = 0; modelIdx < phmmList.count; modelIdx++) {
        numNegativeExamples += (fastaVector.sequence.count / phmmList.phmms[modelIdx].header.modelLength);
    }

    std::cout << "fasta len: " << fastaVector.sequence.count << std::endl;

    const double primaryTpRate = static_cast<double>(numPrimaryTruePositives) / static_cast<double>(numPrimaryPositiveExamples);
    const double complementTpRate = static_cast<double>(numComplementTruePositives) / static_cast<double>(numComplementPositiveExamples);
    const double primaryFpRate = static_cast<double>(numPrimaryFalsePositives) / static_cast<double>(numImpliedSequences - numPrimaryPositiveExamples);
    const double complementFpRate = static_cast<double>(numComplementFalsePositives) / static_cast<double>(numImpliedSequences - numComplementPositiveExamples);

    std::cout << "PTP:\t" << primaryTpRate << "\tPFP:\t" << primaryFpRate << "\tCTP:\t" << complementTpRate <<
        "\tCFP:\t" << complementFpRate << "\n";
    std::cout << "PTP:\t" << numPrimaryTruePositives << "/" << numPrimaryPositiveExamples <<
        "\tPFP:\t" << numPrimaryFalsePositives << "/" << numImpliedSequences - numPrimaryPositiveExamples <<
        "\tCTP:\t" << numComplementTruePositives << "/" << numComplementPositiveExamples <<
        "\tCFP:\t" << numComplementFalsePositives << "/" << numImpliedSequences - numComplementPositiveExamples << std::endl;

    fastaVectorDealloc(&fastaVector);
    p7HmmListDealloc(&phmmList);
}

unique_ptr<vector<vector<PositionListEntry>>> getDnaPositionList(const std::string& positionListFileSrc, const P7HmmList& phmmList, const FastaVector& fastaVector) {
    auto entryList = std::make_unique<vector<vector<PositionListEntry>>>();
    entryList->resize(phmmList.count);

    std::ifstream fileStream(positionListFileSrc);
    std::string line;
    while (std::getline(fileStream, line)) {
        std::istringstream lineStream(line);
        std::string modelLine, modelName, sequenceName, tmp;
        uint64_t positionFrom, positionTo;

        lineStream >> modelLine;
        std::istringstream modelStream(modelLine);
        std::getline(modelStream, modelName, '/');

        lineStream >> sequenceName;
        lineStream >> tmp;
        lineStream >> tmp;
        positionFrom = std::stoi(tmp);
        lineStream >> tmp;
        positionTo = std::stoi(tmp);


        for (uint32_t seqIdx = 0; seqIdx < fastaVector.metadata.count; seqIdx++) {
            size_t seqHeaderLen;
            char* seqHeader;
            fastaVectorGetHeader(&fastaVector, seqIdx, &seqHeader, &seqHeaderLen);

            if (std::string(seqHeader).starts_with(sequenceName)) {

                //find the model idx
                bool foundMatch = false;
                for (uint32_t modelIdx = 0; modelIdx < phmmList.count;modelIdx++) {
                    if (std::string(phmmList.phmms[modelIdx].header.name).ends_with(modelName)) {
                        foundMatch = true;
                        std::string seqNameCopy = sequenceName;
                        (*entryList)[modelIdx].emplace_back(modelName, seqNameCopy, positionFrom, positionTo);
                        break;
                    }
                }
                if (!foundMatch) {
                    std::cout << "unable to find match!" << std::endl;
                }

                break;
            }
        }
    }
    return entryList;
}

std::vector<std::vector<PositionListEntry>> filterDnaPositionLists(const std::vector<std::vector<PositionListEntry>> hmmerHits,
    const std::vector<std::vector<PositionListEntry>> officialPositionList) {

    std::vector<std::vector<PositionListEntry>> filteredEntries;
    filteredEntries.resize(officialPositionList.size());
    for (uint32_t modelIdx = 0; modelIdx < officialPositionList.size();modelIdx++) {
        for (const auto& officialPosition : officialPositionList[modelIdx]) {
            for (const auto& hmmerHit : hmmerHits[modelIdx]) {
                if (officialPosition == hmmerHit) {
                    filteredEntries[modelIdx].push_back(officialPosition);
                    break;
                }
            }
        }
    }
    return filteredEntries;
}


void writeResultsToFile(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complementSeedList,
    const std::string& hmmSrc, const std::string& fastaSrc) {

    FastaVector fastaVector;
    FastaVectorReturnCode fvrc = fastaVectorInit(&fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not init fasta vector in sensitivity check" << std::endl;
        exit(-1);
    }
    fvrc = fastaVectorReadFasta(fastaSrc.c_str(), &fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not read fasta in sensitivity check" << std::endl;
        exit(-1);
    }
    P7HmmList phmmList;
    P7HmmReturnCode phrc = readP7Hmm(hmmSrc.c_str(), &phmmList);
    if (phrc != p7HmmSuccess) {
        std::cerr << "Error: could not read phmm file in sensitivity check" << std::endl;
        exit(-1);
    }

    std::ofstream resultsFileStream;
    resultsFileStream.open("benchmarkResults.txt", std::ios::out);
    resultsFileStream << "#seqName\tmodelName\tseqPosition\tmodelPostiion\tscore\tisComplement\n";


    for (size_t modelIdx = 0; modelIdx < primarySeedList.size(); modelIdx++) {
        //get the model name
        const char* modelName = phmmList.phmms[modelIdx].header.name;
        for (const auto& primaryHit : primarySeedList[modelIdx]) {
            char* sequenceName;
            size_t seqNameLength;
            fastaVectorGetHeader(&fastaVector, primaryHit.sequenceIdx, &sequenceName, &seqNameLength);
            std::string stringSeqName = std::string(sequenceName);
            stringSeqName = stringSeqName.substr(0, stringSeqName.find(' '));
            resultsFileStream << stringSeqName;
            resultsFileStream << "\t";
            resultsFileStream << modelName;
            resultsFileStream << "\t";
            resultsFileStream << primaryHit.sequencePosition;
            resultsFileStream << "\t";
            resultsFileStream << primaryHit.modelPosition;
            resultsFileStream << "\t";
            resultsFileStream << primaryHit.fullScore;
            resultsFileStream << "\t";

            resultsFileStream << "0\n";
        }

        for (const auto& complementHit : complementSeedList[modelIdx]) {
            char* sequenceName;
            size_t seqNameLength;
            fastaVectorGetHeader(&fastaVector, complementHit.sequenceIdx, &sequenceName, &seqNameLength);
            resultsFileStream << sequenceName;
            resultsFileStream << "\t";
            resultsFileStream << modelName;
            resultsFileStream << "\t";
            resultsFileStream << complementHit.sequencePosition;
            resultsFileStream << "\t";
            resultsFileStream << complementHit.modelPosition;
            resultsFileStream << "\t";
            resultsFileStream << complementHit.fullScore;
            resultsFileStream << "\t";
            resultsFileStream << "1\n";
        }
    }

    fastaVectorDealloc(&fastaVector);
    p7HmmListDealloc(&phmmList);
}


void writeHitsToFile(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complementSeedList, const std::string& fastaSrc,
    const std::string hmmSrc, const std::string& outputFileSrc) {

    FastaVector fastaVector;
    FastaVectorReturnCode fvrc = fastaVectorInit(&fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not init fasta vector in sensitivity check" << std::endl;
        exit(-1);
    }
    fvrc = fastaVectorReadFasta(fastaSrc.c_str(), &fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not read fasta in sensitivity check" << std::endl;
        exit(-1);
    }
    P7HmmList phmmList;
    P7HmmReturnCode phrc = readP7Hmm(hmmSrc.c_str(), &phmmList);
    if (phrc != p7HmmSuccess) {
        std::cerr << "Error: could not read phmm file in sensitivity check" << std::endl;
        exit(-1);
    }

    std::ofstream outputStream;
    outputStream.open(outputFileSrc, std::ios::out);

    for (uint32_t modelIdx = 0; modelIdx < primarySeedList.size(); modelIdx++) {
        char* modelName = phmmList.phmms[modelIdx].header.name;
        char* seqName;
        size_t seqNameLen;
        for (const auto& primarySeed : primarySeedList[modelIdx]) {
            fastaVectorGetHeader(&fastaVector, primarySeed.sequenceIdx, &seqName, &seqNameLen);

            std::string stringSeqName = std::string(seqName);
            const auto firstSpaceIdx = stringSeqName.find(" ");
            std::string shortenedSeqName = stringSeqName.substr(0, firstSpaceIdx);

            outputStream << "!!!\t" << shortenedSeqName << "\t" << "TRAIN." << modelName << "\t" <<
                primarySeed.sequencePosition << "\t" << primarySeed.sequencePosition + 10 << "\n";
        }
        for (const auto& complementSeed : complementSeedList[modelIdx]) {
            fastaVectorGetHeader(&fastaVector, complementSeed.sequenceIdx, &seqName, &seqNameLen);
            std::string stringSeqName = std::string(seqName);
            const auto firstSpaceIdx = stringSeqName.find(" ");
            std::string shortenedSeqName = stringSeqName.substr(0, firstSpaceIdx);
            outputStream << "!!!\t" << shortenedSeqName << "\t" << "TRAIN." << modelName << "\t" <<
                complementSeed.sequencePosition + 10 << "\t" << complementSeed.sequencePosition << "\n";
        }
    }
    outputStream.close();


    fastaVectorDealloc(&fastaVector);
    p7HmmListDealloc(&phmmList);
}


struct ScoreEntry {
    ScoreEntry(const float score, const bool isTP, const bool isFP) :
        score(score), isTP(isTP), isFP(isFP) {}
    float score;
    bool isTP;
    bool isFP;
};

void writeAminoScoreList(std::vector<std::vector<NailForge::AlignmentSeed>>& seedList, const std::string& fastaSrc,
    const std::string hmmSrc) {
    FastaVector fastaVector;
    FastaVectorReturnCode fvrc = fastaVectorInit(&fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not init fasta vector in sensitivity check" << std::endl;
        exit(-1);
    }
    fvrc = fastaVectorReadFasta(fastaSrc.c_str(), &fastaVector);
    if (fvrc != FASTA_VECTOR_OK) {
        std::cerr << "Error: could not read fasta in sensitivity check" << std::endl;
        exit(-1);
    }
    P7HmmList phmmList;
    P7HmmReturnCode phrc = readP7Hmm(hmmSrc.c_str(), &phmmList);
    if (phrc != p7HmmSuccess) {
        std::cerr << "Error: could not read phmm file in sensitivity check" << std::endl;
        exit(-1);
    }

    //set a false boolean for each sequence
    std::vector<bool> sequenceHasBeenSeenList;
    sequenceHasBeenSeenList.resize(fastaVector.metadata.count);
    std::fill(sequenceHasBeenSeenList.begin(), sequenceHasBeenSeenList.end(), false);

    std::vector<ScoreEntry> scoreEntryList;

    for (uint32_t modelIdx = 0; modelIdx < seedList.size();modelIdx++) {
        const std::string seedModelName = std::string(phmmList.phmms[modelIdx].header.name);

        for (const auto seed : seedList[modelIdx]) {

            bool weveSeenAHitInThisSequence = sequenceHasBeenSeenList[seed.sequenceIdx];
            char* seqName;
            size_t seqNameLen;
            fastaVectorGetHeader(&fastaVector, seed.sequenceIdx, &seqName, &seqNameLen);
            float score = seed.fullScore;
            bool isTP = false;
            bool isFP = false;
            if (std::string(seqName).starts_with(seedModelName) && !weveSeenAHitInThisSequence) {
                sequenceHasBeenSeenList[seed.sequenceIdx] = true;
                isTP = true;
            }
            else if (std::string(seqName).starts_with("decoy")) {
                isFP = true;
            }
            scoreEntryList.emplace_back(score, isTP, isFP);
        }

    }
    std::sort(scoreEntryList.begin(), scoreEntryList.end(), [](const ScoreEntry& e1, const ScoreEntry& e2) {
        return e1.score > e2.score;
        });

    std::cout << "sorted score data:" << std::endl;
    size_t numTP = 0;
    size_t numFP = 0;
    std::vector<uint32_t> xList;
    std::vector<uint32_t> yList;
    for (const auto& entry : scoreEntryList) {
        if (entry.isTP) {
            numTP++;
        }
        else if (entry.isFP) {
            numFP++;
        }
        yList.push_back(numTP);
        xList.push_back(numFP);
    }

    std::ofstream outStream;
    outStream.open("fig1Data.txt");
    for (uint32_t entryPosition = 0; entryPosition < scoreEntryList.size();entryPosition++) {
        outStream << yList[entryPosition] << "\t" << xList[entryPosition] << std::endl;
    }
    outStream.close();

    fastaVectorDealloc(&fastaVector);
    p7HmmListDealloc(&phmmList);
}