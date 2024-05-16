#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <regex>
#include <cmath>
#include <ctime>
#include <thread>
#include <algorithm> //for std::sort
#include <cstdlib>  // for getenv
#include <unistd.h> // for access
#include <cstring> // for strlen and mmemcpy
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <cstdio>
#include <map>
#include <random>
#include "sam.h"
#include "hts.h"
#include "zlib.h"
#include "edlib.cpp"
#include "igzip_lib.h"
#include "libdeflate.h"
#include "concurrentqueue.h"
#include "report.cpp"

using namespace std;


uint8_t Base[16] = {0,65,67,0,71,0,0,0,84,0,0,0,0,0,0,78};

int  TGSFilter_usage() {
	cout <<""
		"Usage: tgsfilter -i TGS.raw.fq.gz -x ont -o TGS.clean.fq.gz\n"
		" Input/Output options:\n"
		"   -i   <str>   input of bam/fasta/fastq file\n"
		"   -x   <str>   read type (ont|clr|hifi)\n"
		"   -o   <str>   output of fasta/fastq file instead of stdout\n"
		" Basic filter options:\n"
		"   -l   <int>   min length of read to out [1000]\n"
		"   -L   <int>   max length of read to out\n"
		"   -q  <float>  min Phred average quality score\n"
		"   -Q  <float>  max Phred average quality score\n"
		"   -n   <int>   read number for base content check [100000]\n"
		"   -e   <int>   read end length for base content check [150]\n"
		"   -b  <float>  bias (%) of adjacent base content at read end [1]\n"
		"   -5   <int>   trim bases from the 5' end of the read\n"
		"   -3   <int>   trim bases from the 3' end of the read\n"
		" Adapter filter options:\n"
		"   -a   <str>   adapter sequence file \n"
		"   -A           disable reads filter, only for adapter identify\n"
		"   -N   <int>   read number for adapter identify [100000]\n"
		"   -E   <int>   read end length for adapter trim [150]\n"
		"   -m   <int>   min match length for end adapter [15]\n"
		"   -M   <int>   min match length for middle adapter [35]\n"
		"   -T   <int>   extra trim length for middle adpter on both side [50]\n"
		"   -s  <float>  min similarity for end adapter\n"
		"   -S  <float>  min similarity for middle adapter\n"
		"   -D           discard reads with middle adapter instead of split\n"
		" Downsampling options:\n"
		"   -g   <str>   genome size (k/m/g)\n"
		"   -d   <int>   downsample to the desired coverage (requires -g) \n"
		"   -r   <int>   downsample to the desired number of reads \n"
		"   -R  <float>  downsample to the desired fraction of reads \n"
		"   -k   <int>   kmer size for repeat evaluations [11] \n"
		"   -p   <int>   min repeat length of reads [0] \n"
		"   -F           disable reads filter, only for downsampling\n"
		" Other options:\n"
		"   --qc         disable all filter, only for quality control \n"
		"   -f           force FASTA output (discard quality) \n"
		"   -c   <int>   compression level (0-9) for compressed output [6]\n"
		"   -t   <int>   number of threads [16]\n"
		"   -h           show help [v1.11]\n"
		"\n";
	return 1;
}

int compLevel = 6;
int qType =0; //33 or 64 (old) for fastq

class Para_A24 {
	public:
		//
		string InFile;
		string OutFile;
		string TmpOutFile;
		//
		int MinLen;
		int MaxLen;
		float MinQ;
		float MaxQ;
		int BCNum;
		int BCLen;
		float EndBias;
		int HeadTrim;
		int TailTrim;
		//
		string AdapterFile;
		bool ONLYAD;
		int ADNum;
		int EndLen;
		int EndMatchLen;
		int MidMatchLen;
		int ExtraLen;
		float EndSim;
		float MidSim;
		bool discard;
		//
		uint64_t GenomeSize;
		int DesiredDepth;
		int DesiredNum;
		float DesiredFrac;
		bool Downsample;
		bool Filter;
		int Kmer;
		int MinRepeat;
		//
		bool OnlyQC;
		bool FastaOut;
		string readType;
		int n_thread;
		int Infq;
		int Outfq;
		bool OUTGZ;
		
		int ReadLength;

		Para_A24() {
			InFile="";
			OutFile="";
			readType="";
			
			MinLen=1000;
			MaxLen=2147483647;
			MinQ=-1;
			MaxQ=255;
			BCNum=100000;
			BCLen=150;
			EndBias=1;
			HeadTrim=-1;
			TailTrim=-1;
			
			AdapterFile="";
			ONLYAD=false;
			ADNum=100000;
			EndLen=150;
			EndMatchLen=4;
			MidMatchLen=35;
			ExtraLen=50;
			EndSim=0;
			MidSim=0;
			discard=false;

			GenomeSize=0;
			DesiredDepth=0;
			DesiredNum=0;
			DesiredFrac=0;
			Downsample=false;
			Kmer=11;
			MinRepeat=0;
			Filter=true;
			
			OnlyQC=false;
			FastaOut=false;
			n_thread=16;
			Infq=3;
			Outfq=3;
			OUTGZ=false;
			ReadLength=0;
		}
};

inline void  LogLackArg(string &flag) {
	cerr << "Error: Lack Argument for [ -"<<flag<<" ]"<<endl;
}

uint64_t GetGenomeSize(string &genomeSize) {
    std::string numberPart = genomeSize.substr(0, genomeSize.size() - 1);
    char suffix = genomeSize.back();
    double number = std::stod(numberPart);

    uint64_t result = 0;
    if (suffix == 'k' || suffix == 'K') {
        result = static_cast<uint64_t>(number * 1000.0);
    } else if (suffix == 'm' || suffix == 'M') {
        result = static_cast<uint64_t>(number * 1000000.0);
    } else if (suffix == 'g' || suffix == 'G') {
        result = static_cast<uint64_t>(number * 1000000000.0);
    } else {
        cerr << "Error: Genome size should end with k/m/g or K/M/G" << endl;
        return 0; 
    }

    return result;
}

int TGSFilter_cmd(int argc, char **argv, Para_A24 * P2In) {
	if (argc <= 2) {TGSFilter_usage(); return 1;}

	for(int i = 1; i < argc; i++){
		if(argv[i][0] != '-') {
			cerr << "Error: command option error! please check." << endl;
			return 1;
		}

		string flag=argv[i] ;
		flag=regex_replace(flag,"-","");

		//input/output options
		if (flag == "i" ) {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->InFile=argv[i];
		}
		else if (flag == "o" ) {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->OutFile=argv[i];
		}
		else if (flag == "x" ) {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->readType=argv[i];
		}

		//Basic filter options
		else if (flag == "l") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MinLen=atoi(argv[i]);
			if (P2In->MinLen<100){
				P2In->MinLen=100;
			}
		}
		else if (flag == "L") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MaxLen=atoi(argv[i]);
		}
		else if (flag == "q") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MinQ=atof(argv[i]);
		}
		else if (flag == "Q") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MaxQ=atof(argv[i]);
		}
		else if (flag == "n"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->BCNum=atoi(argv[i]);
		}
		else if (flag == "e"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->BCLen=atoi(argv[i]);
		}
		else if (flag == "b"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndBias=atof(argv[i]);
		}
		else if (flag == "5") {
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->HeadTrim=atoi(argv[i]);
		}
		else if (flag == "3"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->TailTrim=atoi(argv[i]);
		}
		
		//Adapter filter options
		else if (flag == "a"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->AdapterFile=argv[i];
		}
		else if (flag == "A"){
			P2In->ONLYAD=true;
		}
		else if (flag == "N"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->ADNum=atoi(argv[i]);
		}
		else if (flag == "E"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndLen=atoi(argv[i]);
		}
		else if (flag == "m"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndMatchLen=atoi(argv[i]);
		}
		else if (flag == "M"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MidMatchLen=atoi(argv[i]);
		}
		else if (flag == "T"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->ExtraLen=atoi(argv[i]);
		}
		else if (flag == "s"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->EndSim=atof(argv[i]);
			if (P2In->EndSim<0.7){
				P2In->EndSim=0.7;
				cerr << "Warning: re set -s to : "<<P2In->EndSim<<endl;
			}
		}
		else if (flag == "S"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MidSim=atof(argv[i]);
			if (P2In->MidSim<0.8){
				P2In->MidSim=0.8;
				cerr << "Warning: reset -S to : "<<P2In->MidSim<<endl;
			}
		}
		else if (flag == "D"){
			P2In->discard=true;
		}

		//Downsampling options
		else if (flag == "g"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			std::string genomeSize = argv[i];
			P2In->GenomeSize=GetGenomeSize(genomeSize);
			if(P2In->GenomeSize == 0){
				return 1;
			}
		}
		else if (flag == "d"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->DesiredDepth=atoi(argv[i]);
		}
		else if (flag == "r"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->DesiredNum=atoi(argv[i]);
		}
		else if (flag == "R"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->DesiredFrac=atof(argv[i]);
		}
		else if (flag == "k"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->Kmer=atoi(argv[i]);
		}
		else if (flag == "p"){
			if(i + 1 == argc) {LogLackArg(flag); return 1;}
			i++;
			P2In->MinRepeat=atoi(argv[i]);
		}
		else if (flag == "F"){
			P2In->Filter=false;
		}

		//Other options
		else if (flag  ==  "qc") {
			P2In->OnlyQC=true;
		}
		else if (flag  ==  "c") {
			if(i + 1 == argc) {LogLackArg(flag) ; return 1;}
			i++;
			compLevel=atoi(argv[i]);
		}
		else if (flag  ==  "f") {
			P2In->FastaOut=true;
		}
		else if (flag  ==  "t") {
			if(i + 1 == argc) {LogLackArg(flag) ; return 1;}
			i++;
			P2In->n_thread=atoi(argv[i]);
		}
		else if (flag  == "help" || flag  == "h") {
			TGSFilter_usage(); return 1;
		}
		else {
			cerr << "Error: UnKnow argument -"<<flag<<endl;
			return 1;
		}
	}

	// check input and output
	if ((P2In->InFile).empty()) {
		cerr<< "Error: lack argument for the must: -i "<<endl;
		exit(-1);
	}else{
		if (access((P2In->InFile).c_str(), 0) != 0) {
			cerr<<"Error: Can't find this file for -i "<<(P2In->InFile)<<endl;
			exit(-1);
		}
	}

	if (P2In->OnlyQC){
		P2In->Filter=false;
	}

	//check read type
	if (P2In->Filter){
		if ((P2In->readType).empty()) {
			cerr<< "Error: lack argument for the must: -x "<<endl;
			exit(-1);
		}else{
			string readType = P2In->readType;
			if (readType == "CLR" || readType == "clr"){
				cerr <<"INFO: read type: PacBio continuous long read (clr)."<<endl;
				P2In->readType="clr";
			}else if (readType == "HIFI" || readType == "hifi"){
				cerr <<"INFO: read type: PacBio highly accurate long reads (hifi)."<<endl;
				P2In->readType = "hifi";
			}else if (readType == "CCS" || readType == "ccs"){
				cerr <<"INFO: read type: PacBio highly accurate long reads (hifi)."<<endl;
				P2In->readType = "hifi";
			}else if (readType == "ONT" || readType == "ont"){
				cerr <<"INFO: read type: NanoPore reads (ont)."<<endl;
				P2In->readType = "ont";
			}else{
				cerr <<"Error: read type should be : clr/hifi/ccs/ont or CLR/HIFI/CCS/ONT."<<endl;
				exit(-1);
			}
		}

		//set read type
		if (P2In->MidSim == 0){
			if (P2In->readType == "hifi"){
				P2In->MidSim = 0.95;
			}else if (P2In->readType == "clr"){
				P2In->MidSim = 0.9;
			}else if (P2In->readType == "ont"){
				P2In->MidSim = 0.9;
			}
		}

		if (P2In->EndSim==0){
			if (P2In->readType == "hifi"){
				P2In->EndSim=0.9;
			}else if (P2In->readType == "clr"){
				P2In->EndSim=0.8;
			}else if (P2In->readType == "ont"){
				P2In->EndSim=0.75;
			}
		}

		cerr <<"INFO: min similarity for middle adapter: "<<P2In->MidSim<<endl;
		cerr <<"INFO: min similarity for end adapter: "<<P2In->EndSim<<endl;
	}
	
	//check downsampling
	if (P2In->DesiredNum > 0 || P2In->DesiredFrac > 0){
		P2In->Downsample=true;
	}else{
		if (P2In->GenomeSize > 0 || P2In->DesiredDepth > 0){
			if (P2In->GenomeSize > 0 && P2In->DesiredDepth > 0){
				P2In->Downsample=true;
			}else if (P2In->GenomeSize > 0 && P2In->DesiredDepth == 0){
				cerr<< "Error: The desired depth was required, along with the genome size!"<<endl;
				exit(-1);
			}else if (P2In->GenomeSize == 0 && P2In->DesiredDepth > 0){
				cerr<< "Error: The genome size was required, along with the desired depth!"<<endl;
				exit(-1);
			}
		}else{
			P2In->Downsample=false;
		}
	}

	if (!(P2In->Filter) && !(P2In->Downsample) && !(P2In->OnlyQC)){
		cerr <<"Error: Please set functional parameters for filter, downsampling or quality control."<<endl;
		exit(-1);
	}
	
	//check threads
	unsigned int maxThreads = std::thread::hardware_concurrency();
	if (maxThreads > 0 && P2In->n_thread > maxThreads - 1){
		P2In->n_thread = maxThreads-1;
		if (P2In->n_thread <= 32){
			cerr <<"Warning: reset -t to: "<<P2In->n_thread<<endl;
		}
	}

	if (P2In->n_thread > 32){
		P2In->n_thread = 32;
		cerr <<"Warning: reset -t to: "<<P2In->n_thread<<endl;
	}
	
	return  0;
	
}

//read fasta or fastq file
#define inBuffSize  1048576
#define outBuffSize 1048576

struct kseq {
        string name;
        string seq;
        string strand;
        string qual;
};

inline bool ends_with(string const & value,  string const & ending) {
	if (ending.size() > value.size()) return false;
	return  equal(ending.rbegin(), ending.rend(), value.rbegin());
}

class FastxReader {
public:
    FastxReader(string filename)
                : mFilename(filename),
                  mFile(nullptr),
                  state(),
	              gz_hdr(),
                  ret(1),
                  remainingData(),
                  combinedData(),
                  isGZ(2),
                  isFastq(2),
                  in_nbytes(0),
                  combinedOffset(0),
                  combinedSize(0),
                  done(false)
    {   
        inBuff = new unsigned char[inBuffSize];
        outBuff = new unsigned char[outBuffSize];
        init();
    }

    ~FastxReader() {
        if (mFile) {
            fclose(mFile);
            mFile = nullptr;
        }
        delete[] inBuff;
        delete[] outBuff;
    }

    kseq* read(){
        if (isFastq == 1){
            return readFastq();
        }else{
            return readFasta();
        }
    }

private:
    void init(){
        mFile = fopen(mFilename.c_str(), "rb");
        if (mFile == nullptr) {
            std::cerr <<"Failed to open file: " << mFilename<<std::endl;
        }

        if (ends_with(mFilename, ".gz")){
            //check gzip header
            isal_gzip_header_init(&gz_hdr);
            isal_inflate_init(&state);
            state.crc_flag = ISAL_GZIP_NO_HDR_VER;
            state.next_in = inBuff;
            state.avail_in = fread(state.next_in, 1, inBuffSize, mFile);
            ret = isal_read_gzip_header(&state, &gz_hdr);
            if (ret != ISAL_DECOMP_OK) {
                cerr << "Error: invalid gzip header found for file: "<<mFilename<<endl;
                exit(-1);
            }

            isGZ = 1;

            if (ends_with(mFilename, ".fastq.gz") || ends_with(mFilename, ".fq.gz")){
                isFastq = 1;
            }else if (ends_with(mFilename, ".fasta.gz") || ends_with(mFilename, ".fa.gz")){
                isFastq = 0;
            }else{
                isFastq = 2;
            }
        }else{
            isGZ = 0;
            if (ends_with(mFilename, ".fastq") || ends_with(mFilename, ".fq")){
                isFastq = 1;
            }else if (ends_with(mFilename, ".fasta") || ends_with(mFilename, ".fa")){
                isFastq = 0;
            }else{
                isFastq = 2;
            }
        }
    }

    void readToBuff(){
        combinedData.clear();
        combinedData.insert(combinedData.end(), remainingData.begin(), remainingData.end());
        
        if (isGZ == 1){
            bool uncompress=false;
            while (!uncompress){
                if(feof(mFile) && state.avail_in==0){
                    done=true;
                    return;
                }

                if (state.block_state != ISAL_BLOCK_FINISH ){

                    if (state.avail_in == 0) {
                        state.next_in = inBuff;
                        state.avail_in = fread(state.next_in, 1, inBuffSize, mFile);
                    }
                    state.next_out = outBuff;
                    state.avail_out = outBuffSize;

                    ret = isal_inflate(&state);
                    if (ret != ISAL_DECOMP_OK) {
                        cerr << "Error: Error encountered while decompressing file: "<<mFilename<<endl;
                        exit(-1);
                    }else{
                        in_nbytes = outBuffSize - state.avail_out;
                        combinedData.insert(combinedData.end(), outBuff, outBuff + (outBuffSize - state.avail_out));
                        remainingData.clear();
                        uncompress=true;
                    }
                } else {
                    if (state.avail_in > 1 && state.next_in[1] != 139){
                        return;
                    }
                    
                    isal_inflate_reset(&state);
                    state.crc_flag = ISAL_GZIP;
                }
            }
        } else if (isGZ == 0){
            if (feof(mFile)){
                done=true;
                return;
            }

            in_nbytes = fread(inBuff, 1, inBuffSize, mFile);
            combinedData.insert(combinedData.end(), inBuff, inBuff + in_nbytes);
            remainingData.clear();
        }

        combinedOffset = 0;
        combinedSize = combinedData.size();

    }

    string getLine(){
        while (true){
            if (combinedOffset < combinedSize){
                const unsigned char *data=combinedData.data() + combinedOffset;
                size_t size = combinedSize - combinedOffset;
                const unsigned char *newline = static_cast<const unsigned char*>(memchr(data, '\n', size));
                
                if (newline != nullptr) {
                    size_t length = newline - data;
                    if (newline > data && newline[-1] == '\r') {
                        length--;
                    }
                    combinedOffset += length + 1;
                    return std::string(reinterpret_cast<const char*>(data), length);
                } else {
                    remainingData.insert(remainingData.end(), combinedData.begin() + combinedOffset, combinedData.end());
                }
            }

            if (!done){
                 readToBuff();
            }else{
                return "";
            }

        }
    }

    kseq* readFastq(){
    
        kseq* ks = new kseq();

        for (int i=0; i<5; i++){
            ks->name= getLine();
            if (!ks->name.empty() && (ks->name[0] == '@')){
                ks->seq = getLine();
                ks->strand = getLine();
                if (ks->strand[0] == '+' && !(ks->seq.empty())){
                    break;
                }
            }
        }

        if(done){
            return nullptr;
        }

        if (ks->name.empty()) {
            cerr <<"Error: input format wrong!"<<endl;
            delete ks;
            return nullptr;
        } else {
            ks->name = ks->name.substr(1);
        }

        ks->qual = getLine();
        if (ks->qual.empty()) {
            cerr << "Error: quality are empty:" << ks->name << endl;
            delete ks;
            return nullptr;
        }

        if (ks->qual.length() != ks->seq.length()) {
            cerr << "warning: sequence and quality have different length:" << ks->name << endl;
            delete ks;
            return nullptr;
        }

        return ks;
    }

    kseq* readFasta(){
        
        kseq* ks = new kseq();

        for (int i=0; i<3; i++){
            ks->name= getLine();
            if (!ks->name.empty() && (ks->name[0] == '>')){
                break;
            }
        }

        if(done){
            return nullptr;
        }

        if (ks->name.empty()) {
            cerr <<"Error: input format wrong!"<<endl;
            delete ks;
            return nullptr;
        } else {
            ks->name = ks->name.substr(1);
        }


        ks->seq = getLine();
        if (ks->seq.empty()) {
            cerr << "Error: sequence are empty:" << ks->name << endl;
            delete ks;
            return nullptr;
        }

        return ks;
    }

    //
    std::string mFilename;
    FILE* mFile;
    struct inflate_state state;
	struct isal_gzip_header gz_hdr;
    unsigned char *inBuff;
    unsigned char *outBuff;
    int ret;

    std::vector<unsigned char> remainingData;
    std::vector<unsigned char> combinedData;

    int isGZ;
    int isFastq;

    size_t in_nbytes;
    size_t combinedOffset;
    size_t combinedSize;
    bool done;

};
//


class DeflateCompress {
public:
    DeflateCompress() {
        mCompressor = libdeflate_alloc_compressor(compLevel);
    }

    ~DeflateCompress() {
        libdeflate_free_compressor(mCompressor);
    }

    bool compressData(const void* input, uint8_t* &out, size_t & outSize) {
        size_t size = strlen(reinterpret_cast<const char*>(input));
        size_t bound = libdeflate_gzip_compress_bound(mCompressor, size);
        out = static_cast<uint8_t*>(malloc(bound));
        outSize = libdeflate_gzip_compress(mCompressor, input, size, out, bound);

        if (outSize == 0) {
            free(out);
            return false;
        } else {
            return true;
        }
    }

private:
    libdeflate_compressor* mCompressor;
};

string GetFileExtension(const string& FilePath) {
	size_t dotPos = FilePath.rfind('.');
	if (dotPos == string::npos) {
		return "";
	} 
	else {
		return FilePath.substr(dotPos + 1);
	}
}

string GetFilePreifx(const string& FilePath){
	string ext = GetFileExtension(FilePath);
	string prefix=FilePath;
	if (ext == "gz") {
		prefix=FilePath.substr(0, FilePath.rfind('.'));
		ext = GetFileExtension(prefix);
	}
	if (ext == "fq" || ext == "fastq" || ext == "fa" || ext == "fasta"){
		prefix=prefix.substr(0, prefix.rfind('.'));
	}else if (ext == "bam" || ext == "sam" || ext == "BAM" || ext == "SAM"){
		prefix=prefix.substr(0, prefix.rfind('.'));
	}
	return prefix;
}

int GetFileType(const string& FilePath) {
	int fqfile=3;
	string ext = GetFileExtension(FilePath);
	if (ext == "gz") {
		string FilePathA = FilePath.substr(0, FilePath.rfind('.'));
		ext = GetFileExtension(FilePathA);
	}

	if ((ext == "fa") || (ext == "fasta")) {
		fqfile = 0;
	} else if ((ext == "fq") || (ext == "fastq")) {
		fqfile = 1;
	}else if ((ext == "sam") || (ext == "SAM") || (ext == "bam") || (ext == "BAM")){
		fqfile = 2;
	}else{
		fqfile = 3;
	}
	return fqfile;
}

char complement[256];
string rev_comp_seq(const string& dna) {
	
	string reverse_complement;
	for (int i = dna.size() - 1; i >= 0; --i) {
		reverse_complement += complement[dna[i]];
	}
	return reverse_complement;
} 

class GetFilterParameterTask {
public:
    GetFilterParameterTask(Para_A24 *P2In, 
						  std::vector<std::string> adapterLib)
                          : P2In(P2In),
						    adapterLib(adapterLib),
                            InPath(P2In->InFile),
                            seqNum(0),
                            maxSeq(0),
                            checkLen(0),
                            minLen(0),
                            seqs5p(),
                            seqs3p(),
                            minQ(255),
							maxQ(0),
							trim5p(0),
							trim3p(0),
							adapter5p(),
							adapter3p(),
							adapterDep5p(0),
							adapterDep3p(0) {}

    int trim5p, trim3p;
    std::string adapter5p, adapter3p;
    float adapterDep5p, adapterDep3p;

    void start() {

        checkLen = P2In->EndLen;
        if (checkLen < P2In->BCLen){
            checkLen = P2In->BCLen;
        }

		if (checkLen < 100){
			checkLen = 100;
		}

		minLen = P2In->MinLen;
		if (minLen < 2 * checkLen){
			minLen = 2 * checkLen;
		}

        maxSeq=P2In->ADNum;
        if (maxSeq<P2In->BCNum){
            maxSeq=P2In->BCNum;
        }

        if ((P2In->Infq)==2){
			read_bam();
		}else if((P2In->Infq)==1 || (P2In->Infq)==0){
			read_fastx();
		}

		if ((P2In->Infq)==1 || (P2In->Infq)==2){
			Get_qType();
		}
		
		if (P2In->Filter){

			std::vector<std::thread> threads;
			
			if ((P2In->HeadTrim)<0){
				threads.emplace_back(&GetFilterParameterTask::CheckBaseContent, this, ref(seqs5p), "5p");
			}
			if ((P2In->TailTrim)<0){
				threads.emplace_back(&GetFilterParameterTask::CheckBaseContent, this, ref(seqs3p), "3p");
			}

			if ((P2In->AdapterFile).empty()){
				threads.emplace_back(&GetFilterParameterTask::adapterSearch, this, ref(seqs5p), "5p");
				threads.emplace_back(&GetFilterParameterTask::adapterSearch, this, ref(seqs3p), "3p");
			}

			for (auto& thread : threads) {
				thread.join();
			}
		}
    }

private:
    int read_fastx() {
		FastxReader reader(InPath);
        kseq* ks = nullptr;
        while ((ks = reader.read()) != nullptr) {
			string seq = ks->seq;
			string qual = ks->qual;
            delete ks;
			int seqLen = seq.length();
            if (seqLen < minLen){
                continue;
            }

            if (seqNum >= maxSeq){
                break;
            }
			seqNum++;

			string seq5p=seq.substr(0, checkLen);
			string seq3p=rev_comp_seq(seq.substr(seqLen-checkLen));
			seqs5p.push_back(seq5p);
			seqs3p.push_back(seq3p);

			qual=qual.substr(0, checkLen);
			for (char q : qual) {
				if(minQ > q) {
					minQ = q;
				}
				if(maxQ < q) {
					maxQ = q;
				}
			}
        }
        return 0;
    }

    int read_bam() {
        htsFile *bamFile = hts_open(InPath.c_str(), "r");
		if (bamFile == nullptr) {
			std::cerr << "Error: Failed to open file: " <<InPath <<endl;
			return 1;
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();

        while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			int seqLen = bamRecord->core.l_qseq;

            if (seqLen < minLen){
                continue;
            }

            if (seqNum >= maxSeq){
                break;
            }
			seqNum++;

            uint8_t *seqChar = bam_get_seq(bamRecord);
			uint8_t *qualChar = bam_get_qual(bamRecord);

            string seq5p, seq3p;

            for (int i = 0; i < checkLen; ++i) {
				char base=Base[bam_seqi(seqChar, i)];
				char qual=qualChar[i] + 33;
				seq5p += base;
				
				if(minQ > qual) {
					minQ = qual;
				}
				if(maxQ < qual) {
					maxQ = qual;
				}
			}

			
			for (int i = (seqLen-checkLen); i < seqLen; ++i) {
				char base=Base[bam_seqi(seqChar, i)];
				seq3p += base;
			}

			seq3p=rev_comp_seq(seq3p);
			seqs5p.push_back(seq5p);
			seqs3p.push_back(seq3p); 
		}
        bam_destroy1(bamRecord);
        bam_hdr_destroy(bamHeader);
        hts_close(bamFile);
        return 0;
    }

    void Get_qType(){
		if(minQ >= 33 &&  minQ <= 78  &&  maxQ >= 33 && maxQ <= 127) {
			qType=33;
		}
		else if (minQ >= 64  &&  minQ <= 108  &&  maxQ >= 64 && maxQ <= 127){
			qType=64;
		}
		else if (minQ < 55) {
			qType=33;
		} else {
			qType=64;
		}

		cerr << "INFO: base quality scoring: Phred"<<qType<<endl;

		minQ-=qType;
		maxQ-=qType;

		if (P2In->MinQ>=0){
			if (P2In->MinQ >= maxQ){
				cerr <<"Warning: max base quality score was: "<<maxQ<<endl;
				cerr <<"INFO: Please reset -q parameter."<<endl;
				exit(-1);
			}
		} else {
			if (maxQ > 10 && P2In->readType=="clr"){
				P2In->MinQ=10;
			}else if (maxQ > 20 && P2In->readType=="hifi"){
				P2In->MinQ=20;
			}else if (maxQ > 10 && P2In->readType=="ont"){
				P2In->MinQ=10;
			}else{
				P2In->MinQ=0;
			}
		}
    }

    void CheckBaseContent(std::vector<std::string> &seqs, string flag){
        std::vector<std::vector<int>> basesNum(checkLen, std::vector<int>(4,0));
        char base;
        for (const std::string& seq : seqs) {
            for (int i = 0; i < seq.length(); i++){
                base=seq[i];
                if (base == 'A' || base=='a'){
                    basesNum[i][0]++;
                }else if (base == 'T' || base=='t'){
                    basesNum[i][1]++;
                }else if (base == 'G' || base=='g'){
                    basesNum[i][2]++;
                }else if (base == 'C' || base=='c'){
                    basesNum[i][3]++;
                }
            }
        }
        
		int maxDiff=(seqNum * P2In->EndBias)/100;
		int lastDis=1;
		
        std::vector<int> canPos;
		int trimLen=0;
        for (int i=1; i < checkLen-1; i++){
			bool leftFlag=false;
			bool rightFlag=false;
			int diff;
			int checkLeftDis=i;
			int checkRightDis=checkLen-i-1;
			if (checkLeftDis > 5){
				checkLeftDis = 5;
			}
			if (checkRightDis > 5){
				checkRightDis = 5;
			}

			for (int j=0; j <basesNum[i].size(); j++){
				for (int x=1; x <= checkLeftDis; x++){
					if (abs(basesNum[i][j] - basesNum[i-x][j]) > maxDiff){
						leftFlag=true;
						break;
					}
				}

				for (int x=1; x <= checkRightDis; x++){
					if (abs(basesNum[i+x][j] - basesNum[i][j]) > maxDiff){
						rightFlag=true;
						break;
					}
				}
			}

			if (leftFlag && rightFlag){
				trimLen=i+1;
			}
        }
        //
		if (trim5p > P2In->BCLen){
			trim5p = P2In->BCLen;
		}

        if (flag=="5p"){
            trim5p=trimLen;
        }else if (flag=="3p"){
            trim3p=trimLen;
        }

    }

    void adapterSearch(std::vector<std::string> &reads, string flag){

		std::unordered_map<int, int> maps;
        float minSim=P2In->MidSim;
		if (minSim < 0.9){
			minSim = 0.9;
		}
		
        for (const auto& ts : reads) {
            int tsLen=ts.length();
			for (int i = 0; i < adapterLib.size(); ++i) {
				string qs=adapterLib[i];
				int qsLen=qs.length();
				int minK=static_cast<int>((1-minSim)*qsLen)+1;

				EdlibAlignResult result_for = edlibAlign(qs.c_str(), qsLen, ts.c_str(), tsLen, 
                            edlibNewAlignConfig(minK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
                if (result_for.status == EDLIB_STATUS_OK){
                    int dist=result_for.editDistance;
                    int length=result_for.alignmentLength;
                    int mlen=length-dist;
                    int numAln=result_for.numLocations;
                    if (numAln>0){
						maps[i]+=mlen;
                    }
                }
                edlibFreeAlignResult(result_for);
			}
		}

		//get max and second index
		std::vector<std::pair<int, int>> vecMap(maps.begin(), maps.end());
		std::sort(vecMap.begin(), vecMap.end(), [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
    		return a.second > b.second;
		});

		string adapter;
		float meanDep=0;
		if (!vecMap.empty()) {
			int maxIndex=vecMap.front().first;
			adapter=adapterLib[maxIndex];
			meanDep = static_cast<float>(maps[maxIndex])/(adapter.length());
		}

        if (flag=="5p"){
			if (meanDep >= 2*minSim){
				adapter5p = adapter;
            	adapterDep5p = meanDep;
			}else{
				adapter5p = "";
            	adapterDep5p = 0;
			}   
        }else if (flag=="3p"){
			if (meanDep >= 2*minSim){
				adapter3p = adapter;
            	adapterDep3p = meanDep;
			}else{
				adapter3p = "";
            	adapterDep3p = 0;
			} 
        }
    }

    Para_A24 *P2In;
	std::vector<std::string> adapterLib;
    std::string InPath, OutPath;
    int seqNum, maxSeq, checkLen, minLen, minQ, maxQ;
    std::vector<std::string> seqs5p, seqs3p;
};

void GetEditDistance(Para_A24 *P2In, const string &query, string &target, 
					int &num5p, int &num3p, int &numMid,
                    std::vector<std::vector<int>> &adapterRegions){
	
	float endSim=P2In->EndSim;
	float midSim=P2In->MidSim;
	
	int qLen=query.length(); // adapter length
	int tLen=target.length(); // reads length

	int end5Len=(P2In->EndLen);
	int end3Len=(P2In->EndLen);

	int extraLen=(P2In->ExtraLen);

	int maxK=qLen-(P2In->MidMatchLen)+1;

	//check middle adapter
	int tsmLen=tLen- end5Len - end3Len;
	if (tsmLen>=qLen){
		string tsm = target.substr(end5Len, tsmLen);
		EdlibAlignResult result = edlibAlign(query.c_str(), qLen, tsm.c_str(), tsmLen, 
						edlibNewAlignConfig(maxK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result.status == EDLIB_STATUS_OK){
			int dist=result.editDistance;
			int numAln=result.numLocations;
			int length=result.alignmentLength;
			int mlen=length - dist;
			if (mlen >= (P2In->MidMatchLen)){
				for (int i=0; i<numAln; i++){
					int ts=result.startLocations[i] + end5Len;
					int te=result.endLocations[i] + end5Len + 1;
					float sim = static_cast<float>(mlen) / qLen;
					
					if (sim >= midSim){
						ts=ts - extraLen;
						te=te + extraLen;
						if (ts<0){ts=0;}
						if (te>tLen){te=tLen;}
						numMid++;
						adapterRegions.push_back({ts, te});
					}
				}
			}
		}
		edlibFreeAlignResult(result);
	}

	//check 5' and 3' end
	int checkLen= (P2In->EndLen) + int(qLen / endSim);
	if (checkLen > tLen){
		checkLen = tLen;
	}
	maxK = qLen-(P2In->EndMatchLen)+1;

	//5' end
	if (checkLen >= 5){
		string ts5 = target.substr(0, checkLen);
		EdlibAlignResult result_ts5 = edlibAlign(query.c_str(), qLen, ts5.c_str(), checkLen, 
					edlibNewAlignConfig(maxK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result_ts5.status == EDLIB_STATUS_OK){
			int dist=result_ts5.editDistance;
			int length=result_ts5.alignmentLength;
			int numAln=result_ts5.numLocations;
			int mlen=length - dist;
			if (mlen >= (P2In->EndMatchLen)){
				for (int i=0; i<numAln; i++){
					int ts=result_ts5.startLocations[i];
					int te=result_ts5.endLocations[i]+1;
					float sim = static_cast<float>(mlen) / qLen;
					if (sim >= endSim){
						num5p++;
						adapterRegions.push_back({0, te});
					}
				}
			}
		}
		edlibFreeAlignResult(result_ts5);
	}

	//3' end
	if (checkLen >= 5){
		string ts3=target.substr(tLen-checkLen);
		EdlibAlignResult result_ts3 = edlibAlign(query.c_str(), qLen, ts3.c_str(), checkLen, 
					edlibNewAlignConfig(maxK, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
		if (result_ts3.status == EDLIB_STATUS_OK){
			int dist=result_ts3.editDistance;
			int length=result_ts3.alignmentLength;
			int numAln=result_ts3.numLocations;
			int mlen=length - dist;
			if (mlen >= (P2In->EndMatchLen)){
				for (int i=0; i<numAln; i++){
					int ts=result_ts3.startLocations[i] + tLen - checkLen;
					int te=result_ts3.endLocations[i]+1 + tLen - checkLen;
					float sim = static_cast<float>(mlen) / qLen;
					if (sim >= endSim){
						num3p++;
						adapterRegions.push_back({ts, tLen});
					}
				}
			}
		}
		edlibFreeAlignResult(result_ts3);
	}
}

std::unordered_set<std::string> adapters;
void adapterMap(Para_A24 * P2In, string &rawSeq, int &rawLen,
				std::vector<std::vector<int>> &keepRegions,
                std::vector<uint64_t> &DropInfo){

	int num5p=0;
	int num3p=0;
	int numMid=0;
	std::vector<std::vector<int>> adapterRegions;

	if (P2In->HeadTrim>0){
		if (P2In->HeadTrim >= rawLen){
			adapterRegions.push_back({0, rawLen});
		}else{
			adapterRegions.push_back({0, P2In->HeadTrim});
		}
	}

	if ((P2In->TailTrim)>0){
		if (P2In->TailTrim >= rawLen){
			adapterRegions.push_back({0, rawLen});
		}else{
			adapterRegions.push_back({rawLen-(P2In->TailTrim), rawLen});
		}
	}
	
	for (const std::string& adapter : adapters) {
		GetEditDistance(P2In, adapter, rawSeq, num5p, num3p, numMid, adapterRegions);
	}

	if (numMid>0 && num5p>0 && num3p>0){
		DropInfo[2]++;
	} else if (numMid>0 && num5p>0 && num3p==0){
		DropInfo[3]++;
	} else if (numMid>0 && num5p==0 && num3p>0){
		DropInfo[4]++;
	} else if (numMid==0 && num5p>0 && num3p>0){
		DropInfo[5]++;
	} else if (numMid>0 && num5p==0 && num3p==0){
		DropInfo[6]++;
	} else if (numMid==0 && num5p>0 && num3p==0){
		DropInfo[7]++;
	} else if (numMid==0 && num5p==0 && num3p>0){
		DropInfo[8]++;
	} else if (numMid==0 && num5p==0 && num3p==0){
		DropInfo[9]++;
	}
	
	if (numMid>0 && P2In->discard){
		DropInfo[10]+=rawLen;
	} else {
		
		std::sort(adapterRegions.begin(), adapterRegions.end(), [](const std::vector<int>& a, const std::vector<int>& b) {
			if (a[0] == b[0]) {
				return a[1] < b[1];
			}
			return a[0] < b[0];
		});

		std::vector<std::vector<int>> mergedRegions;
		for (const auto& region : adapterRegions) {
			if (!mergedRegions.empty() && (mergedRegions.back()[1] >= region[0])) {
				mergedRegions.back()[1] = std::max(mergedRegions.back()[1], region[1]);
			} else {
				mergedRegions.push_back(region);
			}
		}


		int keepLen = 0;
		int dropLen = 0;
		int currentStart = 0;
		if (mergedRegions.size()>=1){
			for (const auto& region : mergedRegions) {
				dropLen = region[1]-region[0];
				DropInfo[10] += dropLen;
				if (dropLen == rawLen){
					DropInfo[11]++;
				}

				if (region[0] > currentStart) {
					keepLen=region[0] - currentStart;
					if (keepLen >= (P2In->MinLen) && keepLen <= (P2In->MaxLen)){
						keepRegions.push_back({currentStart, keepLen});
					}else{
						DropInfo[11]++;
						DropInfo[12]+=keepLen;
					}
				}
				currentStart = region[1];
			}

			if (currentStart < rawLen) {
				keepLen=static_cast<int>(rawLen) - currentStart;
				if (keepLen >= (P2In->MinLen) && keepLen <= (P2In->MaxLen)){
					keepRegions.push_back({currentStart, keepLen});
				}else{
					DropInfo[11]++;
					DropInfo[12]+=keepLen;
				}
			}
		}else{
			if (rawLen >= (P2In->MinLen) && rawLen <= (P2In->MaxLen)){
				keepRegions.push_back({currentStart, rawLen});
			}else{
				DropInfo[11]++;
				DropInfo[12]+=rawLen;
			}
		}
	}
}

double CalcAvgQuality(const string &seq, const string &qual, 
                      std::vector<std::vector<uint64_t>> &baseQual,
                      std::vector<std::vector<uint64_t>> &baseCounts) {
    uint64_t seqLen = seq.length();
    uint64_t qualLen = qual.length();
    if (seqLen == 0 || qualLen == 0 || seqLen != qualLen) {
        return 0.0;
    }

	int vectorSize=int(seqLen/100)+1;
	if (baseQual.size() < vectorSize || baseCounts.size() < vectorSize){
		baseQual.resize(vectorSize, std::vector<uint64_t>(5));
		baseCounts.resize(vectorSize, std::vector<uint64_t>(5));
	}

	uint64_t sumQ = 0;
    int qValue = 0;
	char base;
	int index=0;

    for (uint64_t i = 0; i < seqLen; i++) {
        qValue = qual[i] - qType;
        sumQ += qValue;
		base = seq[i];
		index =int(i/100);

		if (base== 'A' || base=='a'){
			baseCounts[index][0]++;
			baseQual[index][0]+=qValue;
		}else if (base== 'T' || base=='t'){
			baseCounts[index][1]++;
			baseQual[index][1]+=qValue;
		} else if (base== 'G' || base=='g'){
			baseCounts[index][2]++;
			baseQual[index][2]+=qValue;
		} else if (base== 'C' || base=='c'){
			baseCounts[index][3]++;
			baseQual[index][3]+=qValue;
		}
		baseCounts[index][4]++;
		baseQual[index][4]+=qValue;
    }
    return static_cast<double>(sumQ) / seqLen;
}

void Get_5p_base_qual(Para_A24 * P2In, const string &seq, const string &qual, 
                      std::vector<std::vector<uint64_t>> &baseQual,
                      std::vector<std::vector<uint64_t>> &baseCounts) {
    uint64_t seqLen = seq.length();
    uint64_t qualLen = qual.length();
    if (seqLen == 0 || qualLen == 0 || seqLen != qualLen) {
        return;
    }

	int length = P2In->BCLen;
	if (length>seqLen) {
		length=seqLen;
	}

	if (baseQual.size() < length || baseCounts.size() < length){
		baseQual.resize(length, std::vector<uint64_t>(5));
		baseCounts.resize(length, std::vector<uint64_t>(5));
	}

	string headSeq=seq.substr(0,length);
	string headQual=qual.substr(0,length);

	int qValue = 0;
	char base;

    for (uint64_t i = 0; i < length; i++) {
		qValue = headQual[i] - qType;
		base = headSeq[i];

		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
			baseQual[i][0]+=qValue;
		}else if (base== 'T' || base=='t'){
			baseCounts[i][1]++;
			baseQual[i][1]+=qValue;
		} else if (base== 'G' || base=='g'){
			baseCounts[i][2]++;
			baseQual[i][2]+=qValue;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][3]++;
			baseQual[i][3]+=qValue;
		}
		baseCounts[i][4]++;
		baseQual[i][4]+=qValue;
    }
}

void Get_3p_base_qual(Para_A24 * P2In, const string &seq, const string &qual, 
                      std::vector<std::vector<uint64_t>> &baseQual,
                      std::vector<std::vector<uint64_t>> &baseCounts) {
    uint64_t seqLen = seq.length();
    uint64_t qualLen = qual.length();
    if (seqLen == 0 || qualLen == 0 || seqLen != qualLen) {
        return;
    }

	int length = P2In->BCLen;
	if (length>seqLen) {
		length=seqLen;
	}

	if (baseQual.size() < length || baseCounts.size() < length){
		baseQual.resize(length, std::vector<uint64_t>(5));
		baseCounts.resize(length, std::vector<uint64_t>(5));
	}

	string tailSeq=seq.substr(seqLen-length);
	string tailQual=qual.substr(seqLen-length);

	int j = 0;
	int qValue = 0;
	char base;

    for (int i = 0; i < length; i++) {
		j=length-i-1;
		qValue = tailQual[j] - qType;
		base = tailSeq[j];

		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
			baseQual[i][0]+=qValue;
		}else if (base== 'T' || base=='t'){
			baseCounts[i][1]++;
			baseQual[i][1]+=qValue;
		} else if (base== 'G' || base=='g'){
			baseCounts[i][2]++;
			baseQual[i][2]+=qValue;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][3]++;
			baseQual[i][3]+=qValue;
		}
		baseCounts[i][4]++;
		baseQual[i][4]+=qValue;
    }
}

void Get_base_counts(const string &seq, 
                    std::vector<std::vector<uint64_t>> &baseCounts){
    uint64_t seqLen = seq.length();
    if (seqLen == 0){
        return;
    }

	int vectorSize=int(seqLen/100)+1;
	if (baseCounts.size() < vectorSize){
		baseCounts.resize(vectorSize, std::vector<uint64_t>(5));
	}

	char base;
	int index=0;
    for (uint64_t i = 0; i < seqLen; i++){
		base = seq[i];
		index=int(i/100);

		if (base== 'A' || base=='a'){
			baseCounts[index][0]++;
		}else if (base== 'T' || base=='t'){
			baseCounts[index][1]++;
		} else if (base== 'G' || base=='g'){
			baseCounts[index][2]++;
		} else if (base== 'C' || base=='c'){
			baseCounts[index][3]++;
		}
		baseCounts[index][4]++;
    }
}

void Get_5p_base_counts(Para_A24 * P2In, const string &seq, 
                    	std::vector<std::vector<uint64_t>> &baseCounts){
    uint64_t seqLen = seq.length();
    if (seqLen == 0){
        return;
    }

	int length = P2In->BCLen;
	if (length>seqLen) {
		length=seqLen;
	}
	if (baseCounts.size() < length){
		baseCounts.resize(length, std::vector<uint64_t>(5));
	}

	string headSeq=seq.substr(0, length);

	char base;

    for (int i = 0; i < length; i++){
		base = headSeq[i];

		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
		}else if (base== 'T' || base=='t'){
			baseCounts[i][1]++;
		} else if (base== 'G' || base=='G'){
			baseCounts[i][2]++;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][3]++;
		}
		baseCounts[i][4]++;
    }
}

void Get_3p_base_counts(Para_A24 * P2In, const string &seq, 
                    	std::vector<std::vector<uint64_t>> &baseCounts){
    uint64_t seqLen = seq.length();
    if (seqLen == 0){
        return;
    }

	int length = P2In->BCLen;
	if (length>seqLen) {
		length=seqLen;
	}
	if (baseCounts.size() < length){
		baseCounts.resize(length, std::vector<uint64_t>(5));
	}

	string tailSeq=seq.substr(seqLen-length);

	int j = 0;
	char base;

    for (int i = 0; i < length; i++){
		j=length-i-1;
		base = tailSeq[j];

		if (base== 'A' || base=='a'){
			baseCounts[i][0]++;
		}else if (base== 'T' || base=='T'){
			baseCounts[i][1]++;
		} else if (base== 'G' || base=='g'){
			baseCounts[i][2]++;
		} else if (base== 'C' || base=='c'){
			baseCounts[i][3]++;
		}
		baseCounts[i][4]++;
    }
}

std::string newSeqName(const std::string& rawName, int number) {
    std::string newName;
	std::string addNum=":" + std::to_string(number);
	
    bool foundEmpty = false;

    for (char c : rawName) {
        if (std::isspace(c) && !(foundEmpty)) {
			newName += addNum;
			newName += c;
            foundEmpty = true;
		}else{
			newName += c;
		}
    }

    if (!foundEmpty) {
        newName += addNum;
    }

    return newName;
}

int GetKmerCount(const char* seq, int seqLen, int k) {
    std::unordered_set<std::bitset<64>> kmers;

    std::bitset<64> kmer;
    for (int i = 0; i < k; ++i) {
        kmer <<= 2;
        switch (seq[i]) {
            case 'A':
                kmer |= 0;
                break;
            case 'C':
                kmer |= 1;
                break;
            case 'G':
                kmer |= 2;
                break;
            case 'T':
                kmer |= 3;
                break;
            default:
                break;
        }
    }
    kmers.insert(kmer);

    for (int i = k; i < seqLen; ++i) {
        kmer <<= 2;
        switch (seq[i]) {
            case 'A':
                kmer |= 0;
                break;
            case 'C':
                kmer |= 1;
                break;
            case 'G':
                kmer |= 2;
                break;
            case 'T':
                kmer |= 3;
                break;
            default:
                break;
        }
        kmer &= (1ULL << (2 * k)) - 1;
        kmers.insert(kmer);
    }

    int uniqueKmerCount = kmers.size();
    int totalKmerCount = seqLen - k + 1;
    return totalKmerCount - uniqueKmerCount;
}

class TGSFilterTask {
public:
    TGSFilterTask(Para_A24 *P2In)
        : P2In(P2In), 
		  input_queue(), 
		  output_queue(), 
		  outputGz_queue(),
          readNum(0), 
		  filterNum(0), 
		  read_done(false), 
		  filter_done(false),
		  outNum(0), 
		  writeNum(0), 
		  inQueueSize(0), 
		  outQueueSize(0), 
          seqLens(),
		  rawBases(0),
		  cleanBases(0), 
		  rawLens(), 
		  cleanLens(),
		  worker_count(P2In->n_thread),
		  DropInfo(P2In->n_thread, std::vector<uint64_t>(17,0)),
		  rawDiffQualReadsBases(P2In->n_thread, std::vector<uint64_t>(256,0)),
		  cleanDiffQualReadsBases(P2In->n_thread, std::vector<uint64_t>(256,0)),
          rawBaseQual(P2In->n_thread),
		  rawBaseCounts(P2In->n_thread),
		  raw5pBaseQual(P2In->n_thread),
		  raw5pBaseCounts(P2In->n_thread),
		  raw3pBaseQual(P2In->n_thread),
		  raw3pBaseCounts(P2In->n_thread),
	  	  cleanBaseQual(P2In->n_thread),
	  	  cleanBaseCounts(P2In->n_thread),
		  clean5pBaseQual(P2In->n_thread),
		  clean5pBaseCounts(P2In->n_thread),
		  clean3pBaseQual(P2In->n_thread),
		  clean3pBaseCounts(P2In->n_thread) {}

	std::unordered_map<std::string, int> seqLens;
	uint64_t rawBases, cleanBases;
	std::vector<int> rawLens, cleanLens;
	
	std::vector<std::vector<uint64_t>> DropInfo;
	
	std::vector<std::vector<uint64_t>> rawDiffQualReadsBases, cleanDiffQualReadsBases;

	std::vector<std::vector<std::vector<uint64_t>>> rawBaseQual, rawBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> raw5pBaseQual, raw5pBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> raw3pBaseQual, raw3pBaseCounts;

	std::vector<std::vector<std::vector<uint64_t>>> cleanBaseQual, cleanBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> clean5pBaseQual, clean5pBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> clean3pBaseQual, clean3pBaseCounts;

    void start() {

		std::vector<std::thread> read_threads; // read
		if ((P2In->Infq)==2){
			read_threads.emplace_back(&TGSFilterTask::read_bam, this);
		} else if((P2In->Infq)==1 || (P2In->Infq)==0){
			read_threads.emplace_back(&TGSFilterTask::read_fastx, this);
		}

		std::vector<std::thread> filter_threads; // filter
		for (int i = 0; i < worker_count; ++i) {
			filter_threads.emplace_back(&TGSFilterTask::filter_sequence, this, i);
		}

		//
		std::vector<std::thread> output_threads; // output
		if (!(P2In->OnlyQC)){
			if (P2In->OUTGZ && !(P2In->Downsample)){
				output_threads.emplace_back(&TGSFilterTask::write_output_gz, this);
			} else {
				output_threads.emplace_back(&TGSFilterTask::write_output, this);
			}
		}
		
		// Wait for all threads to finish
		for (auto& read_thread : read_threads) {
			read_thread.join();
		}
		for (auto& filter_thread : filter_threads) {
			filter_thread.join();
		}
		for (auto& output_thread : output_threads) {
			output_thread.join();
		}
	}

private:
    int read_fastx() {
		//
		string InPath = P2In->InFile;
		FastxReader reader(InPath);
        kseq* ks = nullptr;
        while ((ks = reader.read()) != nullptr) {
			string name = ks->name;
			string seq = ks->seq;
			string qual = ks->qual;
			int seqLen = seq.length();
			rawBases += seqLen;
			delete ks;

			rawLens.push_back(seqLen);

			input_queue.enqueue(std::make_tuple(name, seq, qual));
			inQueueSize++;
			readNum++;
			if (inQueueSize >= 2 * worker_count){
				this_thread::sleep_for(chrono::milliseconds(1));
			}
            
        }
        read_done = true;
        return 0;
    }

    int read_bam() {
        htsFile *bamFile = hts_open((P2In->InFile).c_str(), "r");
		if (bamFile == nullptr) {
			std::cerr << "Error: Failed to open file: " <<P2In->InFile <<endl;
			return 1;
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();

        while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			string name = bam_get_qname(bamRecord);

			uint8_t *seqChar = bam_get_seq(bamRecord);
			uint8_t *qualChar = bam_get_qual(bamRecord);

			int seqLen = bamRecord->core.l_qseq;
			rawBases += seqLen;
			rawLens.push_back(seqLen);

			string seq;
			string qual;

			for (int i = 0; i < seqLen; ++i) {
				char b = Base[bam_seqi(seqChar, i)];
				char q = qualChar[i] + 33;
				seq += b;
				qual += q;
			}

			input_queue.enqueue(std::make_tuple(name, seq, qual));
			inQueueSize++;
			readNum++;
			if (inQueueSize >= 2 * worker_count){
				this_thread::sleep_for(chrono::milliseconds(1));
			}
		}
        read_done = true;
        bam_destroy1(bamRecord);
        bam_hdr_destroy(bamHeader);
        hts_close(bamFile);
        return 0;
    }

    //
    void filter_sequence(int tid) {
        while (!(read_done && inQueueSize == 0 && readNum == filterNum)) {
			string rawName;
            string rawSeq;
            string rawQual;
			int rawSeqLen = 0;
			int rawQualLen=0;
			
			std::tuple<std::string, std::string, std::string> info;
			if (input_queue.try_dequeue(info)) {
				inQueueSize--;
				rawName = std::get<0>(info);
				rawSeq = std::get<1>(info);
				rawQual = std::get<2>(info);
				rawSeqLen = rawSeq.length();
				rawQualLen = rawQual.length();
			} else {
				this_thread::sleep_for(chrono::milliseconds(1));
			}
			
			if (rawSeqLen > 0){
				if (rawQualLen > 0){
					double rawQuality;
					rawQuality = CalcAvgQuality(rawSeq, rawQual, rawBaseQual[tid], rawBaseCounts[tid]);
					rawDiffQualReadsBases[tid][int(rawQuality)]+=rawSeqLen;
					Get_5p_base_qual(P2In, rawSeq, rawQual, raw5pBaseQual[tid], raw5pBaseCounts[tid]);
					Get_3p_base_qual(P2In, rawSeq, rawQual, raw3pBaseQual[tid], raw3pBaseCounts[tid]);
					if (P2In->Filter){
						if ((rawQuality < (P2In->MinQ)) || (rawQuality > (P2In->MaxQ))){
							DropInfo[tid][0]++;
							DropInfo[tid][1]+=rawSeqLen;
							filterNum++;
							continue;
						}
					}
				} else {
                    Get_base_counts(rawSeq, rawBaseCounts[tid]);
					Get_5p_base_counts(P2In, rawSeq, raw5pBaseCounts[tid]);
					Get_3p_base_counts(P2In, rawSeq, raw3pBaseCounts[tid]);
                }

				std::vector<std::vector<int>> keepRegions;
				if (P2In->Filter){
					adapterMap(P2In, rawSeq, rawSeqLen, keepRegions, DropInfo[tid]);
				}else{
					keepRegions.push_back({0, rawSeqLen});
				}

				string cleanName;
				string cleanSeq;
                string cleanQual;
				string cleanOut;
				int cleanSeqLen;
				int start;
				
				int passNum=1;

				if (keepRegions.size()>0 && !(P2In->OnlyQC)){
					for (const auto& region : keepRegions) {
						start = region[0];
						cleanSeqLen = region[1];
						cleanSeq = rawSeq.substr(start,cleanSeqLen);
						//
						if ((P2In->MinRepeat) > 0){
							int repeatLen = GetKmerCount(cleanSeq.c_str(), cleanSeqLen, P2In->Kmer);
							if (repeatLen < (P2In->MinRepeat)){
								DropInfo[tid][15]++;
								DropInfo[tid][16]+=cleanSeqLen;
								continue;
							}
						}
						//
                        if (rawQualLen > 0){
                            cleanQual=rawQual.substr(start,cleanSeqLen);
							double cleanQuality;
                            cleanQuality = CalcAvgQuality(cleanSeq, cleanQual, cleanBaseQual[tid], cleanBaseCounts[tid]);
							if (P2In->Filter){
								if ((cleanQuality < (P2In->MinQ)) || (cleanQuality > (P2In->MaxQ))){
									DropInfo[tid][13]++;
									DropInfo[tid][14]+=cleanSeqLen;
									continue;
								}
							}
							cleanDiffQualReadsBases[tid][int(cleanQuality)]+=cleanSeqLen;
							Get_5p_base_qual(P2In, cleanSeq, cleanQual, clean5pBaseQual[tid], clean5pBaseCounts[tid]);
							Get_3p_base_qual(P2In, cleanSeq, cleanQual, clean3pBaseQual[tid], clean3pBaseCounts[tid]);
                        }else{
                            Get_base_counts(cleanSeq, cleanBaseCounts[tid]);
							Get_5p_base_counts(P2In, cleanSeq, clean5pBaseCounts[tid]);
							Get_3p_base_counts(P2In, cleanSeq, clean3pBaseCounts[tid]);
                        }

						if (passNum>=2){
							cleanName = newSeqName(rawName, passNum);
						}else{
							cleanName=rawName;
						}
						
						passNum++;
						outNum++;

						if ((P2In->Outfq)==1){
							cleanOut = "@" + cleanName +"\n" + cleanSeq + "\n+\n" + cleanQual + "\n";
							if (P2In->OUTGZ && !(P2In->Downsample)){
								uint8_t *ComData;
								size_t ComSize;
								DeflateCompress GZData;
								if (GZData.compressData(cleanOut.c_str(), ComData, ComSize)) {
									outputGz_queue.enqueue(std::make_tuple(ComData, ComSize, cleanSeqLen));
									outQueueSize++;
								} else {
									free(ComData);
								}
							} else{
								output_queue.enqueue(std::make_tuple(cleanOut, cleanName, cleanSeqLen));
								outQueueSize++;

							}
						} else if((P2In->Outfq)==0){
							cleanOut = ">" + cleanName +"\n" + cleanSeq + "\n";
							if (P2In->OUTGZ && !(P2In->Downsample)){
								uint8_t *ComData;
								size_t ComSize;
								DeflateCompress GZData;
								if (GZData.compressData(cleanOut.c_str(), ComData, ComSize)) {
									outputGz_queue.enqueue(std::make_tuple(ComData, ComSize, cleanSeqLen));
									outQueueSize++;
								} else {
									free(ComData);
								}
							} else {
								output_queue.enqueue(std::make_tuple(cleanOut, cleanName, cleanSeqLen));
								outQueueSize++;
							}
						}

						if (outQueueSize >= 2 * worker_count){
							this_thread::sleep_for(chrono::milliseconds(1));
						}
					}
				}
				filterNum++;
			} 
        }
		filter_done = true;
    }
    //

    void write_output_gz() {
        std::ofstream OUTHGZ((P2In->OutFile).c_str(), std::ios::out | std::ios::binary);
        if (!OUTHGZ.is_open()) {
			std::cerr << "Error: Failed to open file: " <<P2In->OutFile <<endl;
			return ;
        }

        while (!(read_done && filter_done && inQueueSize == 0 && outNum == writeNum && outQueueSize == 0 && readNum == filterNum )) {
			std::tuple<uint8_t*, size_t, int> info;
            if (outputGz_queue.try_dequeue(info)) {
				uint8_t* ComData=std::get<0>(info);
				size_t ComSize=std::get<1>(info);
                int seqLen=std::get<2>(info);
				cleanBases += seqLen;
				cleanLens.push_back(seqLen);
				if (ComData != nullptr && ComSize > 0) {
					OUTHGZ.write(reinterpret_cast<const char*>(ComData), ComSize);
					free(ComData);
					writeNum++;
					outQueueSize--;
				}
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }
        OUTHGZ.close();
    }
    //
    void write_output() {
		if ((P2In->TmpOutFile).empty() && (P2In->OutFile).empty()){
			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
				std::tuple<std::string, std::string, int> info;
				if (output_queue.try_dequeue(info)) {
					string out=std::get<0>(info);
					string name=std::get<1>(info);
					int seqLen=std::get<2>(info);
					cleanBases += seqLen;
					cleanLens.push_back(seqLen);
					seqLens[name] = seqLen;
					std::cout << out;
					writeNum++;
					outQueueSize--;
				} else {
					std::this_thread::sleep_for(std::chrono::milliseconds(1));
				}
			}
		}else{
			string outName;
			if ((P2In->TmpOutFile).empty()){
				outName = (P2In->OutFile);
			}else{
				outName = (P2In->TmpOutFile);
			}

			std::ofstream OUTH(outName.c_str());
			if (!OUTH.is_open()) {
				std::cerr << "Error: Failed to open file: " << outName <<endl;
				return ;
			}

			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
				std::tuple<std::string, std::string, int> info;
				if (output_queue.try_dequeue(info)) {
					string out=std::get<0>(info);
					string name=std::get<1>(info);
					int seqLen=std::get<2>(info);
					cleanBases += seqLen;
					cleanLens.push_back(seqLen);
					seqLens[name] = seqLen;
					OUTH << out;
					writeNum++;
					outQueueSize--; 
				} else {
					std::this_thread::sleep_for(std::chrono::milliseconds(1));
				}
			}
			OUTH.close();
		}
	}

    Para_A24 *P2In;
    moodycamel::ConcurrentQueue<std::tuple<std::string, std::string, std::string>> input_queue;
	moodycamel::ConcurrentQueue<std::tuple<std::string, std::string, uint64_t>> output_queue;
    moodycamel::ConcurrentQueue<std::tuple<uint8_t*, size_t, int>> outputGz_queue;

    std::atomic<bool> read_done;
	std::atomic<bool> filter_done;
	std::atomic<int> readNum;
	std::atomic<int> filterNum;
	std::atomic<int> outNum;
	std::atomic<int> writeNum;
	std::atomic<int> inQueueSize;
	std::atomic<int> outQueueSize;

	int worker_count;
};

class DownSampleTask {
public:
    DownSampleTask(Para_A24 *P2In, string &input, std::unordered_map<std::string, int> &seqLens)
        : P2In(P2In), 
		  InPath(input), 
		  seqLens(seqLens), 
		  worker_count(P2In->n_thread), 
		  totalSize(0),
		  seqNames(), 
		  input_queue(), 
		  output_queue(), 
		  outputGz_queue(),
		  outNum(0), 
		  writeNum(0), 
		  readNum(0), 
		  filterNum(0), 
		  read_done(false), 
		  filter_done(false),
		  inQueueSize(0), 
		  outQueueSize(0),
		  downInNum(0),
		  downInBases(0),
		  downBases(0),
		  downLens(),
		  downDiffQualReadsBases(P2In->n_thread, std::vector<uint64_t>(256,0)),
		  downBaseQual(P2In->n_thread),
	  	  downBaseCounts(P2In->n_thread),
		  down5pBaseQual(P2In->n_thread),
		  down5pBaseCounts(P2In->n_thread),
		  down3pBaseQual(P2In->n_thread),
		  down3pBaseCounts(P2In->n_thread) {}
	
	uint64_t downInNum;
	uint64_t downInBases;
	uint64_t downBases;
	std::vector<int> downLens;
	std::vector<std::vector<uint64_t>> downDiffQualReadsBases;
	std::vector<std::vector<std::vector<uint64_t>>> downBaseQual, downBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> down5pBaseQual, down5pBaseCounts;
	std::vector<std::vector<std::vector<uint64_t>>> down3pBaseQual, down3pBaseCounts;


    void start() {
		//get output names
		if (!(P2In->Filter)){
			if ((P2In->Infq)==2){
				get_bam_SeqLen();
			}else if((P2In->Infq)==1 || (P2In->Infq)==0){
				get_fastx_SeqLen();
			}
		}
		get_reads_name();
		//
		std::vector<std::thread> read_threads; // read
		if (!(P2In->Filter)){
			if ((P2In->Infq)==2){
				read_threads.emplace_back(&DownSampleTask::read_bam, this);
			} else if((P2In->Infq)==1 || (P2In->Infq)==0){
				read_threads.emplace_back(&DownSampleTask::read_fastx, this);
			}
		}else{
			read_threads.emplace_back(&DownSampleTask::read_fastx, this);
		}

		std::vector<std::thread> filter_threads; // filter
		for (int i = 0; i < worker_count; ++i) {
			filter_threads.emplace_back(&DownSampleTask::filter_sequence, this, i);
		}

		//
		std::vector<std::thread> output_threads; // output
		if (P2In->OUTGZ){
			output_threads.emplace_back(&DownSampleTask::write_output_gz, this);
		} else {
			output_threads.emplace_back(&DownSampleTask::write_output, this);
		}

		// Wait for all threads to finish
		
		for (auto& read_thread : read_threads) {
			read_thread.join();
		}
		for (auto& filter_thread : filter_threads) {
			filter_thread.join();
		}
		for (auto& output_thread : output_threads) {
			output_thread.join();
		}
	}


private:
	void get_fastx_SeqLen(){
		FastxReader reader(InPath);
        kseq* ks = nullptr;
        while ((ks = reader.read()) != nullptr) {
			string name = ks->name;
			string seq = ks->seq;
			delete ks;
			int seqLen = seq.length();
			totalSize += seqLen;
			downInNum++;
			downInBases+=seqLen;
			seqLens[name] = seqLen;
		}
	}

	int get_bam_SeqLen(){
		htsFile *bamFile = hts_open(InPath.c_str(), "r");
		if (bamFile == nullptr) {
			cerr <<"Error: Failed to open file: "+InPath<<endl;
			exit(-1);
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();

        while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			string name = bam_get_qname(bamRecord);
			int seqLen = bamRecord->core.l_qseq;
			seqLens[name] = seqLen;
			totalSize += seqLen;
			downInNum++;
			downInBases+=seqLen;
		}
        bam_destroy1(bamRecord);
        bam_hdr_destroy(bamHeader);
        hts_close(bamFile);
		return 0;
	}

	void get_reads_name(){
		std::vector<std::pair<std::string, int>> vec(seqLens.begin(), seqLens.end());
		std::sort(vec.begin(), vec.end(), [](const std::pair<std::string, int>& a, const std::pair<std::string, int>& b) {
    		return a.second > b.second;
		});

		if ((P2In->GenomeSize)>0 && (P2In->DesiredDepth)>0){
			uint64_t addedSize=0;
			uint64_t desiredSize=(P2In->GenomeSize) * (P2In->DesiredDepth);
			for(const auto& pair : vec) {
				seqNames.insert(pair.first);
				addedSize += pair.second;
				downBases += pair.second;
				downLens.push_back(pair.second);
				if (addedSize >= desiredSize){
					break;
				}
			}
		}else if ((P2In->DesiredFrac)>0){
			if (totalSize==0){
				for(const auto& pair : vec) {
					totalSize += pair.second;
				}
			}
			uint64_t addedSize=0;
			uint64_t desiredSize=(P2In->DesiredFrac) * totalSize;
			for(const auto& pair : vec) {
				seqNames.insert(pair.first);
				addedSize += pair.second;
				downBases += pair.second;
				downLens.push_back(pair.second);
				if (addedSize >= desiredSize){
					break;
				}
			}
		}else if ((P2In->DesiredNum)>0){
			int addedNum=0;
			for(const auto& pair : vec) {
				seqNames.insert(pair.first);
				addedNum++;
				downBases += pair.second;
				downLens.push_back(pair.second);
				if (addedNum >= (P2In->DesiredNum)){
					break;
				}
			}
		}
	}

    int read_fastx() {
		FastxReader reader(InPath);
        kseq* ks = nullptr;

        while ((ks = reader.read()) != nullptr) {
			string name = ks->name;
			string seq = ks->seq;
			string qual = ks->qual;
			delete ks;

			if (seqNames.find(name) != seqNames.end()) {
				input_queue.enqueue(std::make_tuple(name, seq, qual));
				inQueueSize++;
				readNum++;
				if (inQueueSize >= 2 * worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			}
        }
		read_done = true;

        return 0;
    }

    int read_bam() {
        htsFile *bamFile = hts_open((P2In->InFile).c_str(), "r");
		if (bamFile == nullptr) {
			std::cerr << "Error: Failed to open file: " <<P2In->InFile <<endl;
			return 1;
		}
		hts_set_opt(bamFile, HTS_OPT_NTHREADS, P2In->n_thread);
		hts_set_log_level(HTS_LOG_OFF);

		bam_hdr_t *bamHeader = sam_hdr_read(bamFile);
		bam1_t *bamRecord = bam_init1();

        while (sam_read1(bamFile, bamHeader, bamRecord) >= 0) {
			string name = bam_get_qname(bamRecord);
			uint8_t *seqChar = bam_get_seq(bamRecord);
			uint8_t *qualChar = bam_get_qual(bamRecord);
			int seqLen = bamRecord->core.l_qseq;

			if (seqNames.find(name) != seqNames.end()) {
				string seq;
				string qual;

				for (int i = 0; i < seqLen; ++i) {
					char b = Base[bam_seqi(seqChar, i)];
					char q = qualChar[i] + 33;
					seq += b;
					qual += q;
				}
				
				input_queue.enqueue(std::make_tuple(name, seq, qual));
				inQueueSize++;
				readNum++;
				if (inQueueSize >= 2 * worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
			}
		}
        read_done = true;
        bam_destroy1(bamRecord);
        bam_hdr_destroy(bamHeader);
        hts_close(bamFile);
        return 0;
    }

    //
    void filter_sequence(int tid) {
        while (!(read_done && inQueueSize == 0 && readNum == filterNum)) {
			string downName;
            string downSeq;
            string downQual;
            int downQualLen=0;
			int downSeqLen = 0;

			std::tuple<std::string, std::string, std::string> info;
			if (input_queue.try_dequeue(info)) {
				inQueueSize--;
				outNum++;
				downName = std::get<0>(info);
				downSeq = std::get<1>(info);
				downQual = std::get<2>(info);
				downSeqLen = downSeq.length();
				downQualLen = downQual.length();
			} else {
				this_thread::sleep_for(chrono::milliseconds(1));
			}
			
			if (downSeqLen > 0){
				if (downQualLen > 0){
					double downQuality;
                    downQuality=CalcAvgQuality(downSeq, downQual, downBaseQual[tid], downBaseCounts[tid]);
					downDiffQualReadsBases[tid][int(downQuality)]+=downSeqLen;
					Get_5p_base_qual(P2In, downSeq, downQual, down5pBaseQual[tid], down5pBaseCounts[tid]);
					Get_3p_base_qual(P2In, downSeq, downQual, down3pBaseQual[tid], down3pBaseCounts[tid]);
                }else{
                    Get_base_counts(downSeq, downBaseCounts[tid]);
					Get_5p_base_counts(P2In, downSeq, down5pBaseCounts[tid]);
					Get_3p_base_counts(P2In, downSeq, down3pBaseCounts[tid]);
                }

				string downOut;
				if ((P2In->Outfq)==1){
					downOut = "@" + downName +"\n" + downSeq + "\n+\n" + downQual + "\n";
					if (P2In->OUTGZ){
						uint8_t *ComData;
						size_t ComSize;
						DeflateCompress GZData;
						if (GZData.compressData(downOut.c_str(), ComData, ComSize)) {
							outputGz_queue.enqueue({ComData, ComSize});
							outQueueSize++;
						} else {
							free(ComData);
						}
					} else{
						output_queue.enqueue(downOut);
						outQueueSize++;

					}
				} else if((P2In->Outfq)==0){
					downOut = ">" + downName +"\n" + downSeq + "\n";
					if (P2In->OUTGZ){
						uint8_t *ComData;
						size_t ComSize;
						DeflateCompress GZData;
						if (GZData.compressData(downOut.c_str(), ComData, ComSize)) {
							outputGz_queue.enqueue({ComData, ComSize});
							outQueueSize++;
						} else {
							free(ComData);
						}
					} else {
						output_queue.enqueue(downOut);
						outQueueSize++;
					}
				}

				if (outQueueSize >= 2 * worker_count){
					this_thread::sleep_for(chrono::milliseconds(1));
				}
				filterNum++;
			}
        }
		filter_done = true;
    }
    //

    void write_output_gz() {
        std::ofstream OUTHGZ((P2In->OutFile).c_str(), std::ios::out | std::ios::binary);
        if (!OUTHGZ.is_open()) {
            throw std::runtime_error("Failed to open output file");
        }

        while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
            std::pair<uint8_t*, size_t> out;
            if (outputGz_queue.try_dequeue(out)) {
                if (out.first != nullptr && out.second > 0) {
                    OUTHGZ.write(reinterpret_cast<const char*>(out.first), out.second);
                    free(out.first);
					writeNum++;
					outQueueSize--;
                }
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
        }

        OUTHGZ.close();
    }
    //
    void write_output() {
		if (!((P2In->OutFile).empty())){
			std::ofstream OUTH((P2In->OutFile).c_str());
        	if (!OUTH.is_open()) {
            	throw std::runtime_error("Failed to open output file");
        	}
			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
            	std::string out;
            	if (output_queue.try_dequeue(out)) {
                	OUTH << out;
					writeNum++;
					outQueueSize--;
            	} else {
                	std::this_thread::sleep_for(std::chrono::milliseconds(1));
            	}
       	 	}

        	OUTH.close();
		}else{
			while (!(read_done && filter_done && readNum == filterNum && outNum == writeNum && outQueueSize == 0 && inQueueSize == 0)) {
            	std::string out;
            	if (output_queue.try_dequeue(out)) {
					std::cout << out;
					writeNum++;
					outQueueSize--;
            	} else {
                	std::this_thread::sleep_for(std::chrono::milliseconds(1));
            	}
        	}
		}
	}

    Para_A24 *P2In;
	string InPath;
	std::unordered_map<std::string, int> seqLens;
	uint64_t totalSize;
	std::unordered_set<std::string> seqNames;
    int worker_count;
    moodycamel::ConcurrentQueue<std::tuple<std::string, std::string, std::string>> input_queue;
    moodycamel::ConcurrentQueue<std::string> output_queue;
    moodycamel::ConcurrentQueue<std::pair<uint8_t*, int>> outputGz_queue;

	std::atomic<bool> read_done;
	std::atomic<bool> filter_done;
	std::atomic<int> readNum;
	std::atomic<int> filterNum;
	std::atomic<int> outNum;
	std::atomic<int> writeNum;
	std::atomic<int> inQueueSize;
	std::atomic<int> outQueueSize;
};


std::string generateRandomString(const int length) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis('a', 'z');

    std::string randomString;
    for (int i = 0; i < length; ++i) {
        randomString += static_cast<char>(dis(gen));
    }

    return randomString;
}

void Get_qual_Dis(const std::vector<std::vector<uint64_t>> &diffQuals,
					uint64_t baseNum, QualDisData &qualDis){
	///merge//
	int maxQual=0;
	std::vector<uint64_t> diffQualsMerge(256,0);
	for (int i = 0; i < diffQuals.size(); i++) {
		for (int j = 0; j < diffQuals[i].size(); j++) {
			uint64_t bases=diffQuals[i][j];
			diffQualsMerge[j]+=bases;
			if (bases>0){
				maxQual=j;
			}
		}
	}
	///ratio///
	qualDis.x.resize(maxQual+1);
	qualDis.y.resize(maxQual+1);
	for (int i = 0; i <= maxQual; i++) {
		qualDis.x[i]=i;
		qualDis.y[i]=static_cast<float>(diffQualsMerge[i]*100)/baseNum;
	}
}

void Get_plot_line_data(Para_A24* P2In,
                        const int &rawMax, 
                        float &rawGC,
                        float &rawMeanQual,
                        const uint64_t &rawBases,
                        LinePlotData &rawReadsQual,
                        LinePlotData &rawBasesContents,
						LinePlotData &raw5pReadsQual,
                        LinePlotData &raw5pBasesContents,
						LinePlotData &raw3pReadsQual,
                        LinePlotData &raw3pBasesContents,
                        const std::vector<std::vector<std::vector<uint64_t>>> &rawBaseQual,
                        const std::vector<std::vector<std::vector<uint64_t>>> &rawBaseCounts,
						const std::vector<std::vector<std::vector<uint64_t>>> &raw5pBaseQual,
						const std::vector<std::vector<std::vector<uint64_t>>> &raw5pBaseCounts,
						const std::vector<std::vector<std::vector<uint64_t>>> &raw3pBaseQual,
                        const std::vector<std::vector<std::vector<uint64_t>>> &raw3pBaseCounts){

    int n_threads=P2In->n_thread;
    uint64_t rawBaseQualSum=0;
    uint64_t rawGCSum=0;
    std::vector<std::string> bases = {"A", "T", "G", "C", "Mean"};
	int baseNum=bases.size();
	//
	uint64_t i=1;
	uint64_t step=100;
	std::vector<uint64_t> indexs;

	while(i < rawMax){
		indexs.push_back(i);
		step=i/20;
		if (step<100){
			step=100;
		}
		
		i+=step;
	}
	indexs.push_back(rawMax);
	uint64_t vecMax=indexs.size()-1;

	std::unordered_map<uint64_t, uint64_t> vecIndexs;
	std::unordered_map<uint64_t, uint64_t> lenIndexs;
    for (uint64_t i = 0; i < vecMax; ++i) {
		uint64_t index=indexs[i];
		vecIndexs[index]=i;
        for (uint64_t x = indexs[i]; x < indexs[i + 1]; ++x) {
            lenIndexs[x] = index;
			//cerr <<x<<" index "<<index<<endl;
        }
    }
	lenIndexs[rawMax]=indexs.back();

	std::vector<std::vector<uint64_t>> rawBaseQualMerge, rawBaseCountsMerge;
	std::vector<std::vector<uint64_t>> raw5pBaseQualMerge, raw5pBaseCountsMerge;
	std::vector<std::vector<uint64_t>> raw3pBaseQualMerge, raw3pBaseCountsMerge;

	int endLen = P2In->BCLen;

	rawBaseCountsMerge.resize(vecMax, std::vector<uint64_t>(baseNum));
	raw5pBaseCountsMerge.resize(endLen, std::vector<uint64_t>(baseNum));
	raw3pBaseCountsMerge.resize(endLen, std::vector<uint64_t>(baseNum));

	rawBaseQualMerge.resize(vecMax, std::vector<uint64_t>(baseNum));
	raw5pBaseQualMerge.resize(endLen, std::vector<uint64_t>(baseNum));
	raw3pBaseQualMerge.resize(endLen, std::vector<uint64_t>(baseNum));

	/////////////////////////////////////////merge data from different threads////////////////////////
	for (int i=0; i<n_threads; i++){
		for (uint64_t j=0; j<rawBaseCounts[i].size(); j++){
			uint64_t lenIndex=lenIndexs[j*100+1];
			uint64_t vecIndex=vecIndexs[lenIndex];
			for (int x=0; x<rawBaseCounts[i][j].size(); x++){
				rawBaseCountsMerge[vecIndex][x]+=rawBaseCounts[i][j][x];
			}
		}
	}

	for (int i=0; i<n_threads; i++){
		for (uint64_t j=0; j < rawBaseQual[i].size();j++){
			uint64_t lenIndex=lenIndexs[j*100+1];
			uint64_t vecIndex=vecIndexs[lenIndex];
			for (int x=0; x < rawBaseQual[i][j].size(); x++){
				rawBaseQualMerge[vecIndex][x]+=rawBaseQual[i][j][x];
			}
		}
	}

	for (int i=0; i<n_threads; i++){
		for (uint64_t j = 0; j < raw5pBaseCounts[i].size(); j++){
			for (int x = 0; x < raw5pBaseCounts[i][j].size(); x++){
				raw5pBaseCountsMerge[j][x] += raw5pBaseCounts[i][j][x];
			}
		}
	}

	for (int i=0; i<n_threads; i++){
		for (uint64_t j = 0; j < raw5pBaseQual[i].size();j++){
			for (int x = 0; x < raw5pBaseQual[i][j].size(); x++){
				raw5pBaseQualMerge[j][x] += raw5pBaseQual[i][j][x];
			}
		}
	}

	for (int i=0; i < n_threads; i++){
		for (uint64_t j = 0; j < raw3pBaseCounts[i].size();j++){
			for (int x = 0; x < raw3pBaseCounts[i][j].size(); x++){
				raw3pBaseCountsMerge[j][x] += raw3pBaseCounts[i][j][x];
			}
		}
	}

	for (int i=0; i<n_threads; i++){
		for (uint64_t j = 0; j < raw3pBaseQual[i].size();j++){
			for (int x = 0; x < raw3pBaseQual[i][j].size(); x++){
				raw3pBaseQualMerge[j][x] += raw3pBaseQual[i][j][x];
				raw3pBaseCountsMerge[j][x] += raw3pBaseCounts[i][j][x];
			}
		}
	}

    ///////////////////////////////////////////////////whole reads/////////////////////////////////////
    rawReadsQual.x.resize(vecMax);
    rawReadsQual.y.resize(baseNum);
    for (int i = 0; i < baseNum; ++i) {
        rawReadsQual.y[i].bases=bases[i];
        rawReadsQual.y[i].data.resize(vecMax);
    }
    
    rawBasesContents.x.resize(vecMax);
    rawBasesContents.y.resize(baseNum-1);
    for (int i = 0; i < baseNum-1; ++i) {
        rawBasesContents.y[i].bases=bases[i];
        rawBasesContents.y[i].data.resize(vecMax);
    }

	uint64_t rawATSum=0;
	for (uint64_t i=0; i < vecMax; i++){
		rawReadsQual.x[i] = indexs[i];
        rawBasesContents.x[i] = indexs[i];
		rawGCSum += rawBaseCountsMerge[i][2];
        rawGCSum += rawBaseCountsMerge[i][3];
		rawATSum += rawBaseCountsMerge[i][0];
		rawATSum += rawBaseCountsMerge[i][1];
		rawBaseQualSum += rawBaseQualMerge[i][baseNum-1];
		for (int j=0; j<baseNum-1; j++){
			if (rawBaseCountsMerge[i][j] > 0){
				rawReadsQual.y[j].data[i] = static_cast<float>(rawBaseQualMerge[i][j])/rawBaseCountsMerge[i][j];
			}else{
				rawReadsQual.y[j].data[i] = 0;
			}

			if (rawBaseCountsMerge[i][baseNum-1] > 0){
				rawBasesContents.y[j].data[i]=static_cast<float>(rawBaseCountsMerge[i][j]*100)/rawBaseCountsMerge[i][baseNum-1];
			}else{
				rawBasesContents.y[j].data[i]=0;
			}
		
		}

		if (rawBaseCountsMerge[i][baseNum-1] >0 ){
			rawReadsQual.y[baseNum-1].data[i]=static_cast<float>(rawBaseQualMerge[i][baseNum-1])/rawBaseCountsMerge[i][baseNum-1];
		}else{
			rawReadsQual.y[baseNum-1].data[i]=0;
		}
	}
	
	rawGC = static_cast<float>(rawGCSum*100)/rawBases; // GC content(%)
    rawMeanQual = static_cast<float>(rawBaseQualSum)/rawBases;

	////////////////////////////////////read 5p ////////////////////////////////////////////////////
	raw5pReadsQual.x.resize(endLen);
    raw5pReadsQual.y.resize(baseNum);
    for (int i = 0; i < baseNum; ++i) {
        raw5pReadsQual.y[i].bases=bases[i];
        raw5pReadsQual.y[i].data.resize(endLen);
    }

    raw5pBasesContents.x.resize(endLen);
    raw5pBasesContents.y.resize(baseNum-1);
    for (int i = 0; i < baseNum-1; ++i) {
        raw5pBasesContents.y[i].bases=bases[i];
        raw5pBasesContents.y[i].data.resize(endLen);
    }

	for (int i=0; i < endLen; i++){
		raw5pReadsQual.x[i] = i+1;
        raw5pBasesContents.x[i] = i+1;
		for (int j=0; j<baseNum-1; j++){
			if (raw5pBaseCountsMerge[i][j] > 0){
				raw5pReadsQual.y[j].data[i]=static_cast<float>(raw5pBaseQualMerge[i][j])/raw5pBaseCountsMerge[i][j];
			}else{
				raw5pReadsQual.y[j].data[i]=0;
			}

			if (raw5pBaseCountsMerge[i][baseNum-1] > 0){
				raw5pBasesContents.y[j].data[i]=static_cast<float>(raw5pBaseCountsMerge[i][j]*100)/raw5pBaseCountsMerge[i][baseNum-1];
			}else{
				raw5pBasesContents.y[j].data[i]=0;
			}    
		}

		if (raw5pBaseCountsMerge[i][baseNum-1] > 0){
			raw5pReadsQual.y[baseNum-1].data[i]=static_cast<float>(raw5pBaseQualMerge[i][baseNum-1])/raw5pBaseCountsMerge[i][baseNum-1];
		}else{
			raw5pReadsQual.y[baseNum-1].data[i]=0;
		}
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////read 3p///////////////////////////////////////////////
	
	raw3pReadsQual.x.resize(endLen);
    raw3pReadsQual.y.resize(baseNum);
    for (int i = 0; i < baseNum; ++i) {
        raw3pReadsQual.y[i].bases=bases[i];
        raw3pReadsQual.y[i].data.resize(endLen);
    }
    
    raw3pBasesContents.x.resize(endLen);
    raw3pBasesContents.y.resize(baseNum-1);
    for (int i = 0; i < baseNum-1; ++i) {
        raw3pBasesContents.y[i].bases=bases[i];
        raw3pBasesContents.y[i].data.resize(endLen);
    }

	for (int i=0; i < endLen; i++){
		raw3pReadsQual.x[i] = i+1;
        raw3pBasesContents.x[i] = i+1;
		for (int j=0; j<baseNum-1; j++){
			if (raw3pBaseCountsMerge[i][j] > 0){
				raw3pReadsQual.y[j].data[i]=static_cast<float>(raw3pBaseQualMerge[i][j])/raw3pBaseCountsMerge[i][j];
			}else{
				raw3pReadsQual.y[j].data[i]=0;
			}
			
			if (raw3pBaseCountsMerge[i][baseNum-1]>0){
				raw3pBasesContents.y[j].data[i]=static_cast<float>(raw3pBaseCountsMerge[i][j]*100)/raw3pBaseCountsMerge[i][baseNum-1];
			}else{
				raw3pBasesContents.y[j].data[i]=0;
			}
		}

		if (raw3pBaseCountsMerge[i][baseNum-1] >0){
			raw3pReadsQual.y[baseNum-1].data[i]=static_cast<float>(raw3pBaseQualMerge[i][baseNum-1])/raw3pBaseCountsMerge[i][baseNum-1];
		}else{
			raw3pReadsQual.y[baseNum-1].data[i]=0;
		}
	}
	
}

void Get_length_Dis(std::vector<int> &rawLens,
						LenDisData &rawLenDis) {
	
	int step=1;
	int i=rawLens.front();
	int rawMax=rawLens.back();
	std::vector<uint64_t> indexs;
	while(i < rawMax){
		indexs.push_back(i);
		step=i/20;
		if (step<100){
			step=100;
		}
		i+=step;
	}
	indexs.push_back(rawMax);

	int vecMax=indexs.size()-1;
	std::unordered_map<int, int> seqIndex;
    for (int i = 0; i < vecMax; ++i) {
        for (int x = indexs[i]; x < indexs[i + 1]; ++x) {
            seqIndex[x] = indexs[i];
        }
    }
	seqIndex[rawMax]=indexs[vecMax-1];

	std::unordered_map<int, int> seqNum;
	for (int len : rawLens) {
		int index=seqIndex[len];
		seqNum[index]++;
	}

	rawLenDis.x.resize(vecMax);
	rawLenDis.y.resize(vecMax);
	for (int i = 0; i < vecMax; ++i) {
		int index=indexs[i];
		rawLenDis.x[i] = index;
		rawLenDis.y[i] = seqNum[index];
		//cerr <<index <<" "<<seqNum[index]<<endl;
	}	
}

int Get_N50(std::vector<int> &rawLens, int &rawNum, uint64_t &rawBases){
	uint64_t rawN50Bases=0;
	for (int i = rawNum - 1; i >= 0; --i) {
		rawN50Bases += rawLens[i];
		if (rawN50Bases >= rawBases/2) {
			return rawLens[i];
		}
	}
	return 0;
}

std::string limitDecimalPlaces(double num, int places) {
    std::string str = std::to_string(num);

    size_t decimalPos = str.find('.');

    if (decimalPos != std::string::npos && str.size() - decimalPos > (places + 1)) {
        str = str.substr(0, decimalPos + places + 1);
    }

    return str;
}

void Get_adapters(Para_A24 * P2In){
	string adapter_name=P2In->AdapterFile;
	FastxReader reader(adapter_name);
	kseq* ks = nullptr;
	while ((ks = reader.read()) != nullptr) {
		string name = ks->name;
		string seq_for = ks->seq;
		string seq_rev=rev_comp_seq(seq_for);
		adapters.insert(seq_for);
		adapters.insert(seq_rev);
		delete ks;
	}

	int num=0;
	for (const std::string& adapter : adapters) {
		num++;
		cerr <<"INFO: input adapter "<<num<<" :"<< adapter << endl;
	}

}

/// ////////////////////////////////////////////////////////////////
int main (int argc, char *argv[ ]) {
	Para_A24 * P2In = new Para_A24;
	int Inflag=0;
	Inflag=TGSFilter_cmd(argc, argv, P2In);
	if(Inflag==1) {
		delete P2In ;
		return 1 ;
	}

	for (int i=0; i<256;i++) {
		complement[i]='N';
	}
	
	complement['A']='T'; complement['G']='C';  
	complement['C']='G'; complement['T']='A';
	complement['a']='t'; complement['g']='c';  
	complement['c']='g'; complement['t']='a';
	complement['M']='K'; complement['R']='Y';  
	complement['W']='W'; complement['S']='S';  
	complement['Y']='R'; complement['K']='M';
	complement['m']='k'; complement['r']='y';  
	complement['w']='w'; complement['s']='s';  
	complement['y']='r'; complement['k']='m';

	std::vector<std::string> adapterLib(22);
	adapterLib[0]="ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT"; // Pacific Biosciences Blunt Adapter
	adapterLib[1]="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT";
	adapterLib[2]="AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA"; // Pacific Biosciences C2 Primer
	adapterLib[3]="TCCTCCTCCTCCGTTAATTTTTTTTTTTTTTTTTT"; 
	adapterLib[4]="AATGTACTTCGTTCAGTTACGTATTGCT"; // Ligation
	adapterLib[5]="AGCAATACGTAACTGAACGAAGTACATT";
	adapterLib[6]="GCAATACGTAACTGAACGAAGT"; // Ligation
	adapterLib[7]="ACTTCGTTCAGTTACGTATTGC";
	adapterLib[8]="GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA"; // Rapid
	adapterLib[9]="TGAAGCGGCGCACGAAAAACGCGAAAGCGTTTCACGATAAATGCGAAAAC";
	adapterLib[10]="GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT"; // 1D^2
	adapterLib[11]="AGGAACCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCC";
	adapterLib[12]="GGAACCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCCAGCAAT"; // 1D^2
	adapterLib[13]="ATTGCTGGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCC";
	adapterLib[14]="TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"; // Ligation(LA) / Native(NA) / Rapid(RA) / Rapid T(RAT) top strand
	adapterLib[15]="AGCAATACGTAACTGAACGAAGTACAGGAAAAAAAA";
	adapterLib[16]="GCAATACGTAACTGAACGAAGTACAGG"; // Ligation Adapter bottom strand
	adapterLib[17]="CCTGTACTTCGTTCAGTTACGTATTGC";
	adapterLib[18]="ACGTAACTGAACGAAGTACAGG"; // Native Adapter bottom strand
	adapterLib[19]="CCTGTACTTCGTTCAGTTACGT";
	adapterLib[20]="CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG"; // cDNA RT Adapter (CRTA)
	adapterLib[21]="CTTGCCTGTCGCTCTATCTTCAGAGGAGAGTCCGCCGCCCGCAAG";

	string InPath=(P2In->InFile);
	P2In->Infq = GetFileType(InPath);
	string prefix=GetFilePreifx(InPath);
	string htmlFileName=prefix + ".html";

	string OutPath;
	if (!(P2In->OutFile).empty() && !(P2In->OnlyQC)) {
		OutPath=(P2In->OutFile);
		htmlFileName = GetFilePreifx(OutPath) + ".html";
		P2In->Outfq = GetFileType(OutPath);
		string ext = GetFileExtension(OutPath);
		if (ext=="gz"){
			P2In->OUTGZ=true;
		}
	}else{
		if (P2In->FastaOut){
			P2In->Outfq = 0;
		}else{
			if (P2In->Infq==2){
				P2In->Outfq = 1;
			}else{
				P2In->Outfq = P2In->Infq;
			}
		}
	}

	if ((P2In->Infq)==3 || (P2In->Outfq)==3) {
		cerr<<"Error: The file name suffix should be '.[fastq|fq|fasta|fa][.gz] or .[sam|bam]'"<<endl;
		if ((P2In->Infq)==3) {
			cerr<<"Error: Please check your input file name: "<<(P2In->InFile)<<endl;
		}
		else if((P2In->Outfq)==3) {
			cerr<<"Error: Please check your output file name: "<<(P2In->OutFile)<<endl;
		}else if ((P2In->Outfq)==2){
			cerr<<"Error: Output file only can be fastq or fasta format: "<<(P2In->OutFile)<<endl;
		}
		return 1;
	}else if ((P2In->Infq)==0 && (P2In->Outfq)==1){
		cerr<<"Error: Fasta format input file can't output fastq format file"<<endl;
		return 1;
	}

	//
	std::unordered_map<std::string, int> seqLens;
	std::vector<int> rawLens;
	std::vector<int> cleanLens;
	std::vector<int> downLens;
	uint64_t rawBases;
	uint64_t cleanBases;
	uint64_t downBases;

	std::vector<std::vector<std::string>> tabInfo(9, std::vector<std::string>(3,"0"));
	LenDisData rawLenDis, cleanLenDis;
	QualDisData rawQualDis, cleanQualDis;

	LinePlotData rawReadsQual, raw5pReadsQual, raw3pReadsQual;
	LinePlotData cleanReadsQual, clean5pReadsQual, clean3pReadsQual;
	
	LinePlotData rawBasesContents, raw5pBasesContents,raw3pBasesContents;
	LinePlotData cleanBasesContents, clean5pBasesContents, clean3pBasesContents;
	
	int rawNum, rawMin, rawMax;
	int cleanNum, cleanMin, cleanMax;

	string downInput;
	GetFilterParameterTask Parameter(P2In, adapterLib);
	Parameter.start();
	
	if (P2In->Filter || P2In->OnlyQC){
		///////////////////////////////////get filter parameter////////////////////////////////////////
		if (P2In->Filter){
			if ((P2In->HeadTrim)<0){
				P2In->HeadTrim=Parameter.trim5p;
			}

			if ((P2In->TailTrim)<0){
				P2In->TailTrim=Parameter.trim3p;
			}
			cerr << "INFO: trim 5' end length: "<<P2In->HeadTrim<<endl;
			cerr << "INFO: trim 3' end length: "<<P2In->TailTrim<<endl;
			cerr << "INFO: min output reads length: "<<P2In->MinLen<<endl;
			if ((P2In->Infq)==1 || (P2In->Infq)==2){
				cerr << "INFO: min Phred average quality score: "<<(P2In->MinQ)<<endl;
			}

			if (!(P2In->AdapterFile).empty()){
				Get_adapters(P2In);
			}else{
				std::string adapter_5p= Parameter.adapter5p;
				std::string adapter_3p= Parameter.adapter3p;
				float depth_5p = Parameter.adapterDep5p;
				float depth_3p = Parameter.adapterDep3p;

				if (depth_5p > 5*depth_3p){
					adapter_3p="";
					depth_3p=0;
				}else if (depth_3p > 5*depth_5p){
					adapter_5p="";
					depth_5p=0;
				}

				cerr <<"INFO: 5' adapter: " << adapter_5p << endl;
				cerr <<"INFO: 3' adapter: " << adapter_3p << endl;

				cerr <<"INFO: mean depth of 5' adapter: "<<depth_5p<<endl;
				cerr <<"INFO: mean depth of 3' adapter: "<<depth_3p<<endl;

				if (P2In->ONLYAD){
					delete P2In ;
					return 0;
				}

				if (!adapter_5p.empty()){
					adapters.insert(adapter_5p);
					adapters.insert(rev_comp_seq(adapter_5p));
				}

				if (!adapter_3p.empty()){
					adapters.insert(adapter_3p);
					adapters.insert(rev_comp_seq(adapter_3p));
				}

				if (adapter_5p.empty() && adapter_3p.empty()){
					if (P2In->readType == "hifi" || P2In->readType == "clr"){
						adapters.insert(adapterLib[0]);
						adapters.insert(adapterLib[1]);
						cerr <<"INFO: set PacBio blunt adapter to trim: "<<adapterLib[0]<<endl;
					}else if (P2In->readType == "ont"){
						adapters.insert(adapterLib[8]);
						adapters.insert(adapterLib[9]);
						cerr <<"INFO: set NanoPore rapid adapter to trim: "<<adapterLib[8]<<endl;
					}
				}
			}

			////////////////////////////////////////////////////////////////////////////////////////////////
			if (P2In->Downsample){
				string rand=generateRandomString(5);
				if ((P2In->Outfq)==0){
					P2In->TmpOutFile = prefix + ".tmp." + rand + ".fa";
				}else if ((P2In->Outfq)==1){
					P2In->TmpOutFile = prefix + ".tmp."+ rand + ".fq";
				}
				downInput=P2In->TmpOutFile;
			}
		}

		TGSFilterTask task(P2In);
    	task.start();
		seqLens=task.seqLens;
		//
		/////////////////////////////////////raw table//////////////////////////////////////////
		
		rawLens=task.rawLens;
		rawBases=task.rawBases;
		rawNum=rawLens.size();
		tabInfo[0][0] = std::to_string(rawNum); // reads number befor filtering
		tabInfo[1][0] = std::to_string(rawBases); // bases number befor filtering
		std::sort(rawLens.begin(), rawLens.end());
		rawMin=rawLens.front();
		rawMax=rawLens.back();
		tabInfo[3][0] = std::to_string(rawMin); // min
		tabInfo[4][0] = std::to_string(rawMax); // max
		tabInfo[5][0] = std::to_string(int(rawBases/rawNum)); // mean
		tabInfo[6][0] = std::to_string(rawLens[int(rawNum/2)]); // median
		tabInfo[7][0] = std::to_string(Get_N50(rawLens, rawNum, rawBases)); // N50
		////////////////////////////////////////////////raw reads base and qual plot/////////////////////////////////////////
		
		float rawGC;
		float rawMeanQual;
		Get_plot_line_data(P2In, rawMax, rawGC, rawMeanQual, rawBases,
                    	rawReadsQual, rawBasesContents,
						raw5pReadsQual, raw5pBasesContents,
						raw3pReadsQual, raw3pBasesContents,
                    	task.rawBaseQual, task.rawBaseCounts,
						task.raw5pBaseQual, task.raw5pBaseCounts,
						task.raw3pBaseQual, task.raw3pBaseCounts);
		tabInfo[2][0] = limitDecimalPlaces(round(rawGC * 1000) / 1000.0, 3); // GC content (% ,string)
		tabInfo[8][0] = limitDecimalPlaces(round(rawMeanQual * 1000) / 1000.0, 3); // mean quality
		Get_length_Dis(rawLens, rawLenDis); // raw Length distribution
		Get_qual_Dis(task.rawDiffQualReadsBases, rawBases, rawQualDis);

		cleanBases=task.cleanBases;
		cleanLens=task.cleanLens;
		cleanNum=cleanLens.size();
		if (!(P2In->OnlyQC) && !(P2In->Downsample)){
			/////////////////////////////////////clean table//////////////////////////////////////////
			tabInfo[0][1] = std::to_string(cleanNum); // reads number after filtering
			tabInfo[1][1] = std::to_string(cleanBases); // bases number after filtering
			std::sort(cleanLens.begin(), cleanLens.end());
			cleanMin=cleanLens.front();
			cleanMax=cleanLens.back();
			tabInfo[3][1] = std::to_string(cleanMin);
			tabInfo[4][1] = std::to_string(cleanMax);
			tabInfo[5][1] = std::to_string(int(cleanBases/cleanNum));
			tabInfo[6][1] = std::to_string(cleanLens[int(cleanNum/2)]);
			tabInfo[7][1] = std::to_string(Get_N50(cleanLens, cleanNum, cleanBases));
			
			////////////////////////////////////////////////clean reads base and qual plot/////////////////////////////////////////
			
			float cleanGC;
			float cleanMeanQual;
			Get_plot_line_data(P2In, cleanMax, cleanGC, cleanMeanQual, cleanBases,
							cleanReadsQual, cleanBasesContents,
							clean5pReadsQual, clean5pBasesContents,
							clean3pReadsQual, clean3pBasesContents,
							task.cleanBaseQual, task.cleanBaseCounts,
							task.clean5pBaseQual, task.clean5pBaseCounts,
							task.clean3pBaseQual, task.clean3pBaseCounts);
			tabInfo[2][1] = limitDecimalPlaces(round(cleanGC * 1000) / 1000.0, 3); // GC content (% ,string)
			tabInfo[8][1] = limitDecimalPlaces(round(cleanMeanQual * 1000) / 1000.0, 3); // mean quality
			Get_length_Dis(cleanLens, cleanLenDis); // clean Length distribution
			Get_qual_Dis(task.cleanDiffQualReadsBases, cleanBases, cleanQualDis);
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::vector<uint64_t> DropInfo(17,0);
		for (int i = 0; i < P2In->n_thread; ++i) {
			for (int j =0; j < 17; ++j){
				DropInfo[j]+=task.DropInfo[i][j];
			}
    	}
		cerr << "INFO: "<< rawNum <<" reads with a total of "<<rawBases<<" bases were input."<<endl;
		if (!(P2In->OnlyQC)){
			cerr << "INFO: "<< DropInfo[0] <<" reads were discarded with "<<DropInfo[1]<<" bases due to low quality."<<endl;
			cerr << "INFO: "<< DropInfo[2] <<" reads have adapter at 5', 3' and middle."<<endl;
			cerr << "INFO: "<< DropInfo[3] <<" reads have adapter at 5' and middle."<<endl;
			cerr << "INFO: "<< DropInfo[4] <<" reads have adapter at 3' and middle."<<endl;
			cerr << "INFO: "<< DropInfo[5] <<" reads have adapter at 5' and 3' end."<<endl;
			cerr << "INFO: "<< DropInfo[6] <<" reads only have adapter at middle."<<endl;
			cerr << "INFO: "<< DropInfo[7] <<" reads only have adapter at 5' end."<<endl;
			cerr << "INFO: "<< DropInfo[8] <<" reads only have adapter at 3' end."<<endl;
			cerr << "INFO: "<< DropInfo[9] <<" reads didn't have any adapter."<<endl;
			cerr << "INFO: "<< DropInfo[10] <<" bases were trimmed due to the adapter or base content bias."<<endl;
			cerr << "INFO: "<< DropInfo[11] <<" reads were discarded with "<<DropInfo[12]<<" bases due to the short length."<<endl;
			cerr << "INFO: "<< DropInfo[13] <<" reads were discarded with "<<DropInfo[14]<<" bases due to low quality after split."<<endl;
			if ((P2In->MinRepeat) > 0){
				cerr << "INFO: "<< DropInfo[15] <<" reads were discarded with "<<DropInfo[16]<<" bases due to short repeat length."<<endl;
			}
			cerr << "INFO: "<< cleanNum <<" reads with a total of "<<cleanBases<<" bases after filtering."<<endl;
			if (!(P2In->Downsample) && !(P2In->OutFile).empty()){
				cerr << "INFO: Filtered reads were written to: "<<P2In->OutFile<<"."<<endl;
			}
		}
	} else {
		downInput=P2In->InFile;
	}

	if (P2In->Downsample){
		DownSampleTask task(P2In, downInput, seqLens);
    	task.start();
		////////////////////////////////////////down table/////////////////////////////////////////
		cleanBases=task.downBases;
		cleanLens=task.downLens;
		cleanNum=cleanLens.size();
		tabInfo[0][1] = std::to_string(cleanNum); // reads number after filtering
		tabInfo[1][1] = std::to_string(cleanBases); // bases number after filtering
		std::sort(cleanLens.begin(), cleanLens.end());
		cleanMin=cleanLens.front();
		cleanMax=cleanLens.back();
		tabInfo[3][1] = std::to_string(cleanMin);
		tabInfo[4][1] = std::to_string(cleanMax);
		tabInfo[5][1] = std::to_string(int(cleanBases/cleanNum));
		tabInfo[6][1] = std::to_string(cleanLens[int(cleanNum/2)]);
		tabInfo[7][1] = std::to_string(Get_N50(cleanLens, cleanNum, cleanBases));
		//
		float cleanGC;
		float cleanMeanQual;
		Get_plot_line_data(P2In, cleanMax, cleanGC, cleanMeanQual, cleanBases,
                    	cleanReadsQual, cleanBasesContents,
						clean5pReadsQual, clean5pBasesContents,
						clean3pReadsQual, clean3pBasesContents,
                    	task.downBaseQual, task.downBaseCounts,
						task.down5pBaseQual, task.down5pBaseCounts,
						task.down3pBaseQual, task.down3pBaseCounts);
		tabInfo[2][1] = limitDecimalPlaces(round(cleanGC * 1000) / 1000.0, 3); // GC content (% ,string)
		tabInfo[8][1] = limitDecimalPlaces(round(cleanMeanQual * 1000) / 1000.0, 2); // mean quality
		Get_length_Dis(cleanLens, cleanLenDis);
		Get_qual_Dis(task.downDiffQualReadsBases, cleanBases, cleanQualDis);
		///////////////////////////////////////////////////////////////////////////////////////
		if (!(P2In->Filter)){
			cerr << "INFO: "<< task.downInNum <<" reads with a total of "<< task.downInBases <<" bases were input."<<endl;
		}
		cerr << "INFO: "<< cleanNum <<" reads with a total of "<<cleanBases<<" bases after downsampling."<<endl;
		if (!(P2In->OutFile).empty()){
			cerr << "INFO: Downsampled reads were written to: "<<P2In->OutFile<<"."<<endl;
		}
	}

	if (!(P2In->TmpOutFile).empty()){
		remove((P2In->TmpOutFile).c_str());
	}

	//get quality control report type
	string qcType;
	if (P2In->Infq==0){
		qcType += "0"; //fasta
	}else if (P2In->Infq==1 || P2In->Infq==2){
		qcType += "1"; //fastq
	}

	if (P2In->OnlyQC){
		qcType += "0";
	}else if (!(P2In->Filter) && P2In->Downsample){
		qcType += "1";
	}else if (P2In->Filter && !(P2In->OnlyQC)){
		qcType += "2";
	}

	PlotData plotData = {
			rawLenDis,
			cleanLenDis,
			rawQualDis,
			cleanQualDis,
			rawReadsQual,
			cleanReadsQual,
			raw5pReadsQual,
			clean5pReadsQual,
			raw3pReadsQual,
			clean3pReadsQual,
			rawBasesContents,
			cleanBasesContents,
			raw5pBasesContents,
			clean5pBasesContents,
			raw3pBasesContents,
			clean3pBasesContents};

	//output html
	ofstream ofs(htmlFileName);
	genHTMLReport(
		ofs, qcType,
		[&tabInfo, &qcType](ofstream &oofs) { genTable(oofs, qcType, tabInfo); },
		[&plotData, &qcType](ofstream &oofs) { genPlotData(oofs, qcType, plotData);}
	);

	ofs.close();
	cerr << "INFO: Quality control report was written to: "<<htmlFileName<<"."<<endl;
	
	delete P2In ;
	return 0;
}