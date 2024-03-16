#include "../src/nailforge.hpp"
#include "../argparse/include/argparse/argparse.hpp"
#include <cctype>
#include <string>

int main(int argc, char** argv) {

    argparse::ArgumentParser parser("benchmark");

    argparse::ArgumentParser createCommand("createfm");
    createCommand.add_description("create the fm index");
    createCommand.add_argument("-f", "--fasta").required().help("specifies the src for the fast file, for generating the awfm index file");
    createCommand.add_argument("-a", "--awFmIndex").required().help("specifies the src for the awfm index file, for either generating or using the index");
    createCommand.add_argument("-r", "suffixArrayRatio").required().help("suffix array compression ratio to use when generating an awfm index.");
    createCommand.add_argument("-A", "--alphabet").required().help("specifies the alphabet for a fasta. values are d for dna, r for rna, a for amino");

    argparse::ArgumentParser searchCommand("search");
    searchCommand.add_argument("-f", "--fasta").required().help("specifies the src for the fast file, for generating the awfm index file");
    searchCommand.add_argument("-a", "--awFmIndex").required().help("specifies the src for the awfm index file, for either generating or using the index");
    searchCommand.add_argument("-h", "--hmmFile").required().help("specifies the src for the hmm file, either for generating a model index or for searching");
    searchCommand.add_argument("-l", "--length").required().help("maximum length a diagonal can accumulate to hit the threshold score in order to generate a hit");
    searchCommand.add_argument("-e", "--extension").required().help("length of the extensions used to validate the hits");
    searchCommand.add_argument("-t", "--thresholdScore").required().help("sets the threshold score in bits for the main diagonal hits.");
    searchCommand.add_argument("-x", "--extensionScore").required().help("the score in bits to hit on the flanking extensions.");
    searchCommand.add_argument("-n", "--numThreads").help("sets the number of threads when openmp is used.");
    searchCommand.arg_argument("-r", "--repeatFilterPartsPerMillion").required().help("float representing maximum number of awfm hits, "
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

    if (parser.present("createfm")) {
        std::string fastaSrc, awfmSrc;
        uint8_t suffixArrayCompressionRatio;
        NailForge::Alphabet alphabet;

        parseCreateOptions(parser, fastaSrc, awfmiSrc, alphabet, suffixArrayCompressionRatio);
        std::cout << "Creating fm index..." << std::endl;
        NailForge::ReturnCode nfrc = NailForge::createFmIndex(fastaSrc.c_str(), awfmSrc.c_str(), alphabet, suffixArrayCompressionRatio);

        if (nfrc != NailForge::ReturnCode::Success) {
            std::cerr << "create fm index returned error code " << nfrc << std::endl;
            exit(-1);
        }
    }
    else if (parser.present("search")) {
        std::string fastaSrc, awfmSrc, hmmSrc;
        NailForge::SearchParams params;
        NailForge::SearchType searchType;
        uint8_t numThreads;
        parseSearchOptions(parser, fastaSrc, awfmSrc, hmmSrc, params, numThreads, searchType);

        std::vector<std::vector<NailForge::AlignmentSeed>> primarySeedList, complimentSeedList;


        std::cout << "filtering with hmm..." << std::endl;
        NailForge::ReturnCode nfrc = NailForge::filterWithHmmFile(hmmFileSrc, fastaSrc, awfmSrc,
            params, searchType, numThreads, primarySeedList, complimentSeedList);

        checkSensitivity(primarySeedList, complimentSeedList, hmmSrc, fastaSrc);

    }
    else {
        std::cout << "no command given, please specify createfm or search" << std::endl;
    }

}


void parseCreateOptions(argparse::ArgumentParser& parser, std::string& fastaSrc, std::string& awfmiSrc, NailForge::Alphabet& alphabet, uint8_t& suffixArrayCompressionRatio) {
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
        std::string alphabetStr = parser.get<std::string>("-A");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get alphabet from parser" << std::endl;
        std::exit(-1);
    }
    if (alphabetStr.size() == 0) {
        switch (std::tolower(alphabetStr[0])) {
        case 'a': alphabet = NailForge::Alphabet::Amino;    break;
        case 'd': alphabet = NailForge::Alphabet::Dna;      break;
        case 'r': alphabet = NailForge::Alphabet::Rna;      break;
        default:
            std::err << "Error, alphabet did not start with d, r, or a." << std::endl;
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
    NailForge::SearchParams& params, NailForge::SearchType& searchType, uint8_t& numThreads) {
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
        params.maximumHitLength = parser.get<uint8_t>("-l");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get main diag length from parser " << std::endl;
        exit(-1);
    }
    try {
        params.flankExtensionLength = parser.get<uint8_t>("-e");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get extension diag length from parser " << std::endl;
        exit(-1);
    }
    try {
        params.extensionThresholdScore = parser.get<float>("-t");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get main diag threshold from parser " << std::endl;
        exit(-1);
    }
    try {
        params.extensionThresholdScore = parser.get<float>("-x");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get extension threshold from parser " << std::endl;
        exit(-1);
    }
    try {
        params.maxSeqHitsPerMillion = parser.get<float>("-r");
    }
    catch (const std::exception& e) {
        std::cerr << "could not parse max seq hits per million from parser" << std::endl;
        exit(-1);
    }
    try {
        numThreads = parser.get<uint8_t>("-n");
    }
    catch (const std::exception& e) {
        std::cerr << "could not get num threads from parser" << std::endl;
        exit(-1);
    }
    try {
        std::string searchTypeString = parser.get<std::string>("-T");
        switch (std::tolower(searchTypeString[0])) {
        case "p": searchType = NailForge::SearchType::Standard;         break;
        case "c": searchType = NailForge::SearchType::DualStrand;       break;
        case "b": searchType = NailForge::SearchType::ComplimentStrand; break;
        case 'd': searchType = NailForge::SearchType::ComplimentStrand; break;
        default:
            std::cerr << "search type was not primary, compliment, or dual-strand" << std::endl;
            exit(-1)
        }
    }
    catch (const std::exception& e) {
        std::cerr << "could not parse search type option from parser" << std::endl;
        exit(-1);
    }
}

void checkSensitivity(const std::vector<std::vector<NailForge::AlignmentSeed>>& primarySeedList,
    const std::vector<std::vector<NailForge::AlignmentSeed>>& complimentSeedLiset,
    const std::string& hmmFileSrc, const std::string& fastaFileSrc) {


}