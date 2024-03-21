#include "nailforge.hpp"
#include "PhmmProcessor/PhmmProcessor.hpp"
#include "StringTree/StringTree.hpp"
#include <string>
#include <iostream>

namespace NailForge {
    bool allModelsInHmmListSameAlphabet(P7HmmList* phmmList);

    NailForge::ReturnCode createFmIndex(const char* fastaFileSrc, const char* fmIndexFileSrc,
        const NailForge::Alphabet& alphabet, const uint8_t suffixArrayCompressionRatio) {

        AwFmIndex* fmIndex;
        struct AwFmIndexConfiguration awfmConfig;
        awfmConfig.suffixArrayCompressionRatio = suffixArrayCompressionRatio;
        awfmConfig.keepSuffixArrayInMemory = true;
        awfmConfig.storeOriginalSequence = false;



        if (alphabet == Alphabet::Dna) {
            awfmConfig.kmerLengthInSeedTable = 12;  //16*4^12 = 268MB
            awfmConfig.alphabetType = AwFmAlphabetDna;
        }
        else if (alphabet == Alphabet::Rna) {
            awfmConfig.kmerLengthInSeedTable = 12;  //16*4^12 = 268MB
            awfmConfig.alphabetType = AwFmAlphabetRna;
        }
        else {
            awfmConfig.kmerLengthInSeedTable = 5;   //16*20^6 = 51MB
            awfmConfig.alphabetType = AwFmAlphabetAmino;
        }

        enum AwFmReturnCode awfmRc = awFmCreateIndexFromFasta(&fmIndex,
            &awfmConfig, fastaFileSrc, fmIndexFileSrc);

        switch (awfmRc) {
        case AwFmFileWriteOkay:                 return ReturnCode::Success;
        case AwFmSuccess:                       return ReturnCode::Success;
        case AwFmNullPtrError:                  return ReturnCode::GeneralFailure;
        case AwFmFileOpenFail:                  return ReturnCode::FileNotFound;
        case AwFmAllocationFailure:             return ReturnCode::MemoryError;
        case AwFmFileAlreadyExists:             return ReturnCode::FileAlreadyExists;
        case AwFmSuffixArrayCreationFailure:    return ReturnCode::GeneralFailure;
        case AwFmFileWriteFail:                 return ReturnCode::FileAlreadyExists;
        default:                                return ReturnCode::GeneralFailure;
        }
    }


    NailForge::ReturnCode filterWithHmmFile(const char* hmmFileSrc, const char* fastaFileSrc, const char* fmIndexFileSrc,
        const NailForge::SearchParams& params, const SearchType searchType, const uint8_t numThreads,
        std::vector<std::vector<AlignmentSeed>>& primarySeedList, std::vector<std::vector<AlignmentSeed>>& complimentSeedList) {

        omp_set_num_threads(numThreads);

        //load the fm index

        AwFmIndex* fmIndex;
        FastaVector fastaVector;
        P7HmmList phmmList;
        enum AwFmReturnCode awfmRc = awFmReadIndexFromFile(&fmIndex, fmIndexFileSrc, true);//always keep the SA in memory
        switch (awfmRc) {
        case AwFmFileAlreadyExists:
            std::cerr << "awfm index file already exists, could not write." << std::endl;
            return NailForge::ReturnCode::GeneralFailure;
        case AwFmFileFormatError:
            std::cerr << "AwFmIndex file src does not seem to point to an AwFmIndex file." << std::endl;
            return NailForge::ReturnCode::FileFormatError;
        case AwFmAllocationFailure:
            std::cerr << "failed to allocate memory for the fm index" << std::endl;
            return NailForge::ReturnCode::AllocationFailure;
        }


        //load the fasta
        fastaVectorInit(&fastaVector);
        FastaVectorReturnCode fvrc = fastaVectorReadFasta(fastaFileSrc, &fastaVector);
        switch (fvrc) {
        case FASTA_VECTOR_ALLOCATION_FAIL:
            awFmDeallocIndex(fmIndex);
            std::cerr << "could not allocate memory for the sequence" << std::endl;
            return NailForge::ReturnCode::AllocationFailure;
        case FASTA_VECTOR_FILE_OPEN_FAIL:
            awFmDeallocIndex(fmIndex);
            std::cerr << "failed to  read fasta file." << std::endl;
            return NailForge::ReturnCode::FileReadError;
        }


        //load the phmmList
        P7HmmReturnCode p7hmmRc = readP7Hmm(hmmFileSrc, &phmmList);
        switch (p7hmmRc) {
        case p7HmmFormatError:
            fastaVectorDealloc(&fastaVector);
            awFmDeallocIndex(fmIndex);
            std::cerr << "file src for .hmm file does not seem to be a hmm file" << std::endl;
            return ReturnCode::FileFormatError;
        case p7HmmAllocationFailure:
            fastaVectorDealloc(&fastaVector);
            awFmDeallocIndex(fmIndex);
            std::cerr << "Failed to allocate memory for the hmm list" << std::endl;
            return ReturnCode::AllocationFailure;
        case p7HmmFileNotFound:
            fastaVectorDealloc(&fastaVector);
            awFmDeallocIndex(fmIndex);
            std::cerr << ".hmm file not found" << std::endl;
            return ReturnCode::FileNotFound;
        }

        if (!allModelsInHmmListSameAlphabet(&phmmList)) {
            std::cerr << "error: all models in a given .hmm file must be the same alphabet." << std::endl;
            fastaVectorDealloc(&fastaVector);
            awFmDeallocIndex(fmIndex);
            p7HmmListDealloc(&phmmList);
            return ReturnCode::AlphabetMismatch;
        }
        //this seems like it'll never happen, but better to have a check for it.
        if (phmmList.count == 0) {
            std::cerr << "error: .hmm file must contain at least one model" << std::endl;
            fastaVectorDealloc(&fastaVector);
            awFmDeallocIndex(fmIndex);
            p7HmmListDealloc(&phmmList);
            return ReturnCode::GeneralFailure;
        }

        NailForge::Alphabet alphabet;
        switch (phmmList.phmms[0].header.alphabet) {
        case P7HmmReaderAlphabetAmino:  alphabet = NailForge::Alphabet::Amino;  break;
        case P7HmmReaderAlphabetDna:    alphabet = NailForge::Alphabet::Dna;    break;
        case P7HmmReaderAlphabetRna:    alphabet = NailForge::Alphabet::Rna;    break;
        default:    std::cerr << "error: hmm file alphabet not supported. only Amino, Rna, and Dna are supported." << std::endl;
            fastaVectorDealloc(&fastaVector);
            awFmDeallocIndex(fmIndex);
            p7HmmListDealloc(&phmmList);
            return ReturnCode::AlphabetMismatch;
        }

        if (!alphabetsMatch(phmmList, fmIndex)) {
            std::cerr << "one or more phmms in the given .hmm file didn't match the alphabet of the fm index." << std::endl;
            fastaVectorDealloc(&fastaVector);
            awFmDeallocIndex(fmIndex);
            p7HmmListDealloc(&phmmList);
            return ReturnCode::AlphabetMismatch;
        }

        primarySeedList.resize(phmmList.count);
        complimentSeedList.resize(phmmList.count);

        uint64_t modelsCompleted = 0;
        // std::cout << "model " << 0 << " / " << phmmList.count << std::endl;
#pragma omp parallel for
        for (uint32_t modelIdx = 0; modelIdx < phmmList.count; modelIdx++) {
            const auto& phmm = phmmList.phmms[modelIdx];
            std::vector<float> matchScores = NailForge::PhmmProcessor::toFloatMatchScores(phmm);

            if (searchType != NailForge::SearchType::ComplimentStrand) {
                NailForge::StringTree::findSeeds(fmIndex, fastaVector, phmm, modelIdx,
                    matchScores, primarySeedList[modelIdx], params, false);
            }
            if (searchType != NailForge::SearchType::Standard) {
                NailForge::StringTree::findSeeds(fmIndex, fastaVector, phmm, modelIdx,
                    matchScores, complimentSeedList[modelIdx], params, true);

            }
            // std::cout << "model "<< modelIdx<< " / "<<phmmList.count<< std::endl;
        }

        fastaVectorDealloc(&fastaVector);
        awFmDeallocIndex(fmIndex);
        p7HmmListDealloc(&phmmList);
        return NailForge::ReturnCode::Success;
    }




    std::string_view returnCodeDescription(const NailForge::ReturnCode rc) {
        switch (rc) {
        case ReturnCode::Success:           return "success";
        case ReturnCode::AllocationFailure: return "allocation failure";
        case ReturnCode::AlphabetMismatch:  return "alphabet mismatch";
        case ReturnCode::FileAlreadyExists: return "file already exists";
        case ReturnCode::FileFormatError:   return "file format error";
        case ReturnCode::FileNotFound:      return "file not found";
        case ReturnCode::FileReadError:     return "file read error";
        case ReturnCode::FileWriteError:    return "file write error";
        case ReturnCode::FmIndexError:      return "fm index error";
        case ReturnCode::GeneralFailure:    return "general failure";
        case ReturnCode::MemoryError:       return "memory exception";
        case ReturnCode::NotImplemented:    return "functionality not implemented";
        default:                            return "no description given for this return code";
        }
    }

    bool allModelsInHmmListSameAlphabet(P7HmmList* phmmList) {
        if (phmmList->count == 0) {
            return true;
        }
        P7Alphabet firstModelAlphabet = phmmList->phmms[0].header.alphabet;

        for (uint32_t i = 1; i < phmmList->count; i++) {
            if (phmmList->phmms[i].header.alphabet != firstModelAlphabet) {
                return false;
            }
        }
        return true;
    }


}
