/* The MIT License

   Copyright (c) 2023- by Huiyang Yu, Weiming He, Chunmei Shi.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#ifndef FQ_KmerSplit_H_
#define FQ_KmerSplit_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
#include <thread>
#include <algorithm>
#include <malloc.h>
#include <atomic>
#include <unordered_map>
#include <map>
#include <cstdio>
#include <vector>
#include <zlib.h>
#include "./gzstream.c"
#include "./kseq.h"
#include "./comm.hpp"
#include "./M1_work.hpp"
#include "./kc-c4-window.c"

typedef  long long LLongA;

using namespace std;
//KSEQ_INIT(gzFile, gzread)

int  print_usage_FqSplit()
{
	cout <<""
		"Usage: hfkreads -1 PE1.fq.gz -2 PE2.fq.gz -o OutFrefix\n"
		" Input/Output options:\n"
		"   -1	<str>   paired-end fasta/q file1\n"
		"   -2	<str>   paired-end fasta/q file2\n"
		"   -s	<str>   single-end fasta/q\n"
		"   -o	<str>   prefix of output file\n"
		" Filter options:\n"
		"   -b	<int>   min base quality [0]\n"
		"   -q	<int>   min average base quality [20]\n"
		"   -l	<int>   min length of read [half]\n"
		"   -r	<float> max unknown base (N) ratio [0.1]\n"
		"   -k	<int>   kmer length [31]\n"
		"   -w	<int>   window size [5]\n"
		"   -m	<int>   min kmer count for high freq kmer [3] \n"
		"   -x	<int>   min count of read with high freq kmer [5]\n"
		"   -n	<int>   read number to use [1000000]\n"
		"   -a	        use all the read number\n"
		" Other options:\n"
		"   -c           compress the outPut File\n"
		"   -f           outPut the KmerFre File\n"
		"   -A           keep output quality info\n"
		"   -v           PE use SE read to remove PCR\n"
		"   -t           thread to run [4]\n"
		"   -h           show help [v2.02]\n"
		"\n";
	return 1;
}
//
int parse_cmd_FqSplit(int argc, char **argv , Para_A24 * P2In  )
{
	if (argc <=3  ) {print_usage_FqSplit();return 0;}

	int err_flag = 0;

	for(int i = 1; i < argc || err_flag; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "Error: command option error! please check." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");
		//Input/Output options
		if (flag  == "1" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InFq1=argv[i];
		}
		else if (flag  == "2")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InFq2=argv[i];
		}
		else if (flag  ==  "s" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->InSeFq=argv[i];
		}
		else if (flag  ==  "o" )
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->OutFq1=argv[i];
		}
		//Filter low quality
		else if (flag  ==  "b")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->MinBaseQ=atoi(argv[i]);
		}
		else if (flag  ==  "q")
		{
			if(i + 1 == argc) { LogLackArg(flag) ; return 0;}
			i++;
			P2In->AverQ=atoi(argv[i]);
		}
		else if (flag  ==  "l")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->HalfReadLength=atoi(argv[i]);
			if  ( (P2In->HalfReadLength)<11 )
			{
				P2In->HalfReadLength=11;
				cerr<<"Warings: -l should >= 11, we set it to 11\n";
			}
		}
		else if (flag  ==  "r")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->N_Ration=atof(argv[i]);
		}
		//Filter low freq reads
		else if (flag  ==  "k")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->Kmer=atoi(argv[i]);
			if ( (P2In->Kmer) >31 || (P2In->Kmer) <  9 )
			{
				cerr<<"Warning -k  [9,31] ,we modify it to be  31\n ";
				P2In->Kmer=31;
			}
		}
		else if (flag  ==  "w")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->Windows=atoi(argv[i]);
			if  ( (P2In->Windows)<1 )
			{
				P2In->Windows=1;
				cerr<<"Warnings: -w must be >=1, we change it to 1 \n";
			}
		}
		else if (flag  ==  "m")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->MinCount=atoi(argv[i]);
			if  ((P2In->MinCount)<1)
			{
				cerr<<"Warnings: -m should >= 1, we set it to 1\n";
				(P2In->MinCount)=1;
			}
		}
		else if (flag  =="x")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->MinReadKmerCount=atoi(argv[i]);
			if( (P2In->MinReadKmerCount)<1 )  { P2In->MinReadKmerCount=1 ;}
		}
		else if (flag  ==  "n")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			P2In->ReadNumber=atoi(argv[i]);
		}
		else if (flag  == "A")
		{
			P2In->OutFa=false;
		}
		//other options
		else if (flag  == "v")
		{
			P2In->PCRA=true;
		}
		else if (flag  == "c")
		{
			P2In->OUTGZ=true;
		}
		else if (flag  == "f")
		{
			P2In->KmerStatOut=true;
		}
		else if (flag  == "a")
		{
			P2In->ReadNumber=LONG_MAX;
			//P2In->allRead=true;
		}
		else if (flag  ==  "t")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			n_thread=atoi(argv[i]);
		}
		else if (flag  == "help" || flag  == "h")
		{
			print_usage_FqSplit();return 0;
		}
		//
		else if (flag  ==  "u")
		{
			if(i + 1 == argc) { LogLackArg( flag ) ; return 0;}
			i++;
			VECMAX=atoi(argv[i]);
		}
		else
		{
			cerr << "Error: UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}

	if ( (P2In->InSeFq).empty() && (P2In->InFq1).empty() )
	{
		cerr<< "Error: -1/-2 or -s lack argument for the must"<<endl;
		return 0;
	}

	if ( (P2In->OutFq1).empty() )
	{
		cerr<< "Error: -o lack argument for the must"<<endl;
		return 0;
	}

	if  ( (P2In->InFq1).empty())
	{
		return 2 ;
	}
	else if ( (P2In->InFq2).empty() || (P2In->InFq1).empty() )
	{
		cerr<< "Error: -1 and -1 must set together"<<endl;
		return 0;
	}
	else
	{
		return  1 ;
	}
}

int RunFQFilterPE (Para_A24 * P2In,  vector<std::string>  & FilePath)
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAQQQ ;
	vector <string>  AAAIII ;


	vector <string>  BBBSSS ;
	vector <string>  BBBQQQ ;
	vector <string>  BBBIII ;


	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE > (P2In->ReadNumber)) 
	{
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
		}

	AAASSS.resize(BATCH_SIZE+2);
	AAAQQQ.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	BBBSSS.resize(BATCH_SIZE+2);
	BBBQQQ.resize(BATCH_SIZE+2);
	BBBIII.resize(BATCH_SIZE+2);

	int A=0;

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	bool  *PASS =new bool [BATCH_SIZE];

	igzstream INA ((P2In->InFq1).c_str(),ifstream::in);
	igzstream INB ((P2In->InFq2).c_str(),ifstream::in);
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	INB.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	string ID_1 ,seq_1,temp_1,Quly_1 ;
	string ID_2 ,seq_2,temp_2,Quly_2 ;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+"_APE_tmp.fa";
	string outputFileBPE=OUT+"_BPE_tmp.fa";
	if (!P2In->OutFa)
	{
		outputFileAPE=OUT+"_APE_tmp.fq";
		outputFileBPE=OUT+"_BPE_tmp.fq";
	}
	FilePath.push_back(outputFileAPE);
	FilePath.push_back(outputFileBPE);

	OUTIO OUTHanDle ;
	OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
	OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
	int CountFQ=0;

	if (P2In->OutFa)
	{

		for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);

			getline(INB,ID_2);
			getline(INB,seq_2);
			getline(INB,temp_2);
			getline(INB,Quly_2);

			if (ID_1.empty()) {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;

			BBBSSS[CountFQ]=seq_2;
			BBBQQQ[CountFQ]=Quly_2;
			BBBIII[CountFQ]=ID_2;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterFQPE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{						
						AAAIII[j][0]='>';
						BBBIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					}
				}
				CountFQ=0;
			}
		}


		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }

				threads.push_back(std::thread(FilterFQPE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
			}
			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{
					AAAIII[j][0]='>';
					BBBIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
				}
			}
			CountFQ=0;
		}
	}
	else
	{
		for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);

			getline(INB,ID_2);
			getline(INB,seq_2);
			getline(INB,temp_2);
			getline(INB,Quly_2);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;

			BBBSSS[CountFQ]=seq_2;
			BBBQQQ[CountFQ]=Quly_2;
			BBBIII[CountFQ]=ID_2;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterFQPE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }

				threads.push_back(std::thread(FilterFQPE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ),std::ref(BBBSSS),std::ref(BBBQQQ)));
			}
			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
				}
			}
			CountFQ=0;
		}

	}

	OUTHanDle.OUTAPE.close();
	OUTHanDle.OUTBPE.close();

	if (!INA.eof()) {  getline(INA,ID_1); if (INA.eof()) { A++;} }
	if (INA.eof()) {A--;	cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;	}

	INA.close();
	INB.close();
	delete [] Start ;
	delete [] End;
	delete [] PASS;


	return 0;
}

int RunFQFilterSE (Para_A24 * P2In,  vector<std::string>  & FilePath)
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAQQQ ;
	vector <string>  AAAIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if  (BATCH_SIZE  > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAQQQ.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);


	int A=0;

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];


	bool  *PASS =new bool [BATCH_SIZE];

	igzstream INA ((P2In->InFq1).c_str(),ifstream::in);
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	string ID_1 ,seq_1,temp_1,Quly_1 ;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+"_SE_tmp.fa";
	if (!P2In->OutFa)
	{
		outputFileAPE=OUT+"_SE_tmp.fq";
	}
	FilePath.push_back(outputFileAPE);

	OUTIO OUTHanDle ;
	OUTHanDle.OUTAPE.open(outputFileAPE.c_str());

	int CountFQ=0;

	if (P2In->OutFa)
	{

		for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);


			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;


			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{						
						AAAIII[j][0]='>';
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
			}
			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRSE(P2In,PASS,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{
					AAAIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
				}
			}
			CountFQ=0;
		}

	}
	else
	{

		for (A=0 ; A<(P2In->ReadNumber) && (!INA.eof()) ; A++)
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRSE(P2In,PASS,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS[j])
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterFQSE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(AAAQQQ)));
			}
			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRSE(P2In,PASS,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
				}
			}
			CountFQ=0;
		}
	}

	OUTHanDle.OUTAPE.close();

	if (!INA.eof()) {  getline(INA,ID_1); if (INA.eof()) { A++;} }
	if (INA.eof()) {A--;	cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;	}

	INA.close();
	delete [] Start ;
	delete [] End;
	delete [] PASS;
	return 0;
}

int RunFAFilterPE (Para_A24 * P2In,  vector<std::string>  & FilePath)
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAIII ;

	vector <string>  BBBSSS ;
	vector <string>  BBBIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if  (BATCH_SIZE  > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}
	AAASSS.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	BBBSSS.resize(BATCH_SIZE+2);
	BBBIII.resize(BATCH_SIZE+2);

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool  *PASS =new bool [BATCH_SIZE];

	gzFile fpAA;
	kseq_t *seqAA;

	gzFile fpBB;
	kseq_t *seqBB;

	fpAA = gzopen((P2In->InFq1).c_str(), "r");
	seqAA = kseq_init(fpAA);

	fpBB = gzopen((P2In->InFq2).c_str(), "r");
	seqBB = kseq_init(fpBB);

	int AA ; int A=0;
	string ID_1,ID_2;
	string seqStrAA,seqStrBB;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+"_APE_tmp.fa";
	string outputFileBPE=OUT+"_BPE_tmp.fa";

	FilePath.push_back(outputFileAPE);
	FilePath.push_back(outputFileBPE);

	OUTIO OUTHanDle ;
	OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
	OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
	int CountFQ=0;

	for (A=0 ; A<(P2In->ReadNumber) && ( (AA = kseq_read(seqAA)) >= 0)  ; A++)
	{
		ID_1=(seqAA->name.s);
		seqStrAA=(seqAA->seq.s);

		AAASSS[CountFQ]=seqStrAA;
		AAAIII[CountFQ]=ID_1;
		AA = kseq_read(seqBB);

		ID_2=(seqBB->name.s);
		seqStrBB=(seqBB->seq.s);
		BBBSSS[CountFQ]=seqStrBB;
		BBBIII[CountFQ]=ID_2;
		CountFQ++;

		if (CountFQ==BATCH_SIZE)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				threads.push_back(std::thread(FilterFAPE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}

			threads.clear();
			RmPCRPE(P2In,PASS,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{						
					AAAIII[j][0]='>';
					BBBIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
				}
			}
			CountFQ=0;
		}
	}

	if (CountFQ!=0)
	{
		for (int i = 0; i < n_thread; i++)
		{
			Start[i]=i*BinWind;
			End[i]=Start[i]+BinWind;
			if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
			threads.push_back(std::thread(FilterFAPE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
		}

		for (auto& thread : threads)
		{
			thread.join();
		}
		threads.clear();
		RmPCRPE(P2In,PASS,CountFQ,AAASSS,BBBSSS);

		for (int j = 0; j < CountFQ; j++)
		{
			if (PASS[j])
			{
				AAAIII[j][0]='>';
				BBBIII[j][0]='>';
				OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
				OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
			}
		}
		CountFQ=0;
	}

	OUTHanDle.OUTAPE.close();
	OUTHanDle.OUTBPE.close();

	if  (A==(P2In->ReadNumber) )
	{
		if  ((AA = kseq_read(seqAA)) >= 0) {  }
		else
		{
			cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
		}
	}
	else
	{
		cout<<"INFO: ALL reads "<<A*2<<" are read done"<<endl;
	}

	kseq_destroy(seqAA);
	gzclose(fpAA);

	kseq_destroy(seqBB);
	gzclose(fpBB);

	delete [] Start ;
	delete [] End;
	delete [] PASS;

	return 0;
}

int RunFAFilterSE (Para_A24 * P2In,  vector<std::string>  & FilePath)
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	bool  *PASS =new bool [BATCH_SIZE];

	gzFile fpAA;
	kseq_t *seqAA;

	fpAA = gzopen((P2In->InFq1).c_str(), "r");
	seqAA = kseq_init(fpAA);

	int AA ; int A=0;
	string ID_1;
	string seqStrAA;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+"_SE_tmp.fa";

	FilePath.push_back(outputFileAPE);

	OUTIO OUTHanDle ;
	OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
	int CountFQ=0;

	for (A=0 ; A<(P2In->ReadNumber) && ( (AA = kseq_read(seqAA)) >= 0)  ; A++)
	{

		ID_1=(seqAA->name.s);
		seqStrAA=(seqAA->seq.s);

		AAASSS[CountFQ]=seqStrAA;
		AAAIII[CountFQ]=ID_1;

		CountFQ++;

		if (CountFQ==BATCH_SIZE)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				threads.push_back(std::thread(FilterFASE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}

			threads.clear();
			RmPCRSE(P2In,PASS,CountFQ,AAASSS);


			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS[j])
				{						
					AAAIII[j][0]='>';
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
				}
			}
			CountFQ=0;
		}
	}

	if (CountFQ!=0)
	{
		for (int i = 0; i < n_thread; i++)
		{
			Start[i]=i*BinWind;
			End[i]=Start[i]+BinWind;
			if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
			threads.push_back(std::thread(FilterFASE,P2In,PASS,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
		}

		for (auto& thread : threads)
		{
			thread.join();
		}
		threads.clear();
		RmPCRSE(P2In,PASS,CountFQ,AAASSS);

		for (int j = 0; j < CountFQ; j++)
		{
			if (PASS[j])
			{
				AAAIII[j][0]='>';
				OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
			}
		}
		CountFQ=0;
	}

	OUTHanDle.OUTAPE.close();

	if  (A==(P2In->ReadNumber) )
	{
		if  ((AA = kseq_read(seqAA)) >= 0) {  }
		else
		{
			cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
		}
	}
	else
	{
		cout<<"INFO: ALL reads "<<A<<" are read done"<<endl;
	}

	kseq_destroy(seqAA);
	gzclose(fpAA);

	delete [] Start ;
	delete [] End;
	delete [] PASS;

	return 0;
}

void GetMinCount(Para_A24 * P2In,  const kc_c4x_t *h  )
{

	if(((P2In->MinCount)==0) || (P2In->KmerStatOut))
	{
		hist_aux_t a;
		uint64_t cnt[256];
		int i, j;
		a.h = h;

		CALLOC(a.cnt, n_thread);
		kt_for(n_thread, worker_hist, &a, 1<<h->p);
		for (i = 0; i < 256; ++i) { cnt[i] = 0;}
		for (j = 0; j < n_thread; ++j)
		{
			for (i = 0; i < 256; ++i)
			{
				cnt[i] += a.cnt[j].c[i];
			}
		}
		free(a.cnt);

		uint64_t  Sum=0;
		for (i = 1; i < 256; ++i)
		{
			Sum+=cnt[i]*i;
		}
		uint64_t  Half=Sum/2;

		if ((P2In->KmerStatOut))
		{
			string KmerFre=(P2In->OutFq1)+".KmerFre.gz";
			ogzstream OUTStat (KmerFre.c_str());
			for (i = 1; i < 256; ++i)
			{
				OUTStat<<i<<"\t"<<cnt[i]<<"\n";
			}
			OUTStat.close();
		}
		else
		{
			Sum=0;
			for (i = 1; i < 256 &&  Sum<Half; ++i)
			{
				Sum+=cnt[i]*i;
			}

			if ((P2In->MinCount)==0) { (P2In->MinCount)=i;}
		}
	}

	if ((P2In->MinCount)<2) {P2In->MinCount=2;}
	cout<<"INFO: min kmer count is set to :"<<(P2In->MinCount)<<endl;
}

int  FilterLowHitPE(Para_A24 * P2In,  const kc_c4x_t *h , bool  * PASS1, bool  * PASS2,  int & Start,int &   End ,vector <string> & AAASSS,vector <string> & BBBSSS)
{
	bool  QPASS_AA;
	bool  QPASS_BB;
	int CountAA=0;
	int CountBB=0;
	//return ;
	for (int ii=Start; ii<End; ii++)
	{

		CountAA=ReadHitNum( h , P2In->Kmer,P2In->Windows,P2In->MinCount,AAASSS[ii]);
		CountBB=ReadHitNum( h , P2In->Kmer,P2In->Windows,P2In->MinCount,BBBSSS[ii]);
		
		if (CountAA<(P2In->MinReadKmerCount))
		{
			PASS1[ii]=false;
		}
		else
		{
			PASS1[ii]= true;
		}

		if (CountBB<(P2In->MinReadKmerCount))
		{

			PASS2[ii]=false;
		}
		else
		{
			PASS2[ii]= true;
		}

	}

	return 1;
}

void FilterLowHitSE(Para_A24 * P2In, const kc_c4x_t *h , bool  * PASS1,  int & Start,int &   End ,vector <string> & AAASSS)
{
	bool  QPASS_AA;
	int CountAA=0;

	for (int ii=Start; ii<End; ii++)
	{

		CountAA=ReadHitNum( h , P2In->Kmer,P2In->Windows,P2In->MinCount,AAASSS[ii]);
		if (CountAA<(P2In->MinReadKmerCount))
		{

			PASS1[ii]=false;
		}
		else
		{
			PASS1[ii]= true;
		}

	}


}

int RunFQ2FQ_PEOUT ( Para_A24 * P2In,  vector<std::string>  & FilePath, const kc_c4x_t *h )
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAQQQ ;
	vector <string>  AAAIII ;

	vector <string>  BBBSSS ;
	vector <string>  BBBQQQ ;
	vector <string>  BBBIII ;


	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if  (BATCH_SIZE  > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAQQQ.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	BBBSSS.resize(BATCH_SIZE+2);
	BBBQQQ.resize(BATCH_SIZE+2);
	BBBIII.resize(BATCH_SIZE+2);


	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];


	bool  *PASS1 =new bool [BATCH_SIZE];
	bool  *PASS2 =new bool [BATCH_SIZE];




	igzstream INA (FilePath[0].c_str(),ifstream::in); // ifstream  + gz
	igzstream INB (FilePath[1].c_str(),ifstream::in); // ifstream  + gz
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	INB.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);


	string ID_1 ,seq_1,temp_1,Quly_1;
	string ID_2 ,seq_2,temp_2,Quly_2;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+"_pe_1.fq";
	string outputFileBPE=OUT+"_pe_2.fq";
	string outputFileASE=OUT+"_se_1.fq";
	string outputFileBSE=OUT+"_se_2.fq";

	if (P2In->OUTGZ)
	{

		OUTIOGZ OUTHanDle ;
		outputFileAPE=OUT+"_pe_1.fq.gz";
		outputFileBPE=OUT+"_pe_2.fq.gz";
		outputFileASE=OUT+"_se_1.fq.gz";
		outputFileBSE=OUT+"_se_2.fq.gz";

		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;

		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		int CountFQ=0;

		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);

			getline(INB,ID_2);
			getline(INB,seq_2);
			getline(INB,temp_2);
			getline(INB,Quly_2);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;

			BBBSSS[CountFQ]=seq_2;
			BBBQQQ[CountFQ]=Quly_2;
			BBBIII[CountFQ]=ID_2;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{

				if (PASS1[j]  & PASS2[j] )
				{						
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
					OUTHanDle.PE++;
				}
				else if (PASS1[j])
				{
					OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					OUTHanDle.SEAA++;
				}
				else if (PASS2[j])
				{
					OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
					OUTHanDle.SEBB++;
				}
			}
			CountFQ=0;
		}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();
	}
	else
	{
		OUTIO OUTHanDle ;

		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;
		int CountFQ=0;

		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);

			getline(INB,ID_2);
			getline(INB,seq_2);
			getline(INB,temp_2);
			getline(INB,Quly_2);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;

			BBBSSS[CountFQ]=seq_2;
			BBBQQQ[CountFQ]=Quly_2;
			BBBIII[CountFQ]=ID_2;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j]  & PASS2[j] )
				{						
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
					OUTHanDle.PE++;
				}
				else if (PASS1[j])
				{
					OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					OUTHanDle.SEAA++;
				}
				else if (PASS2[j])
				{
					OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n+\n"<<BBBQQQ[j]<<"\n";
					OUTHanDle.SEBB++;
				}
			}
			CountFQ=0;
		}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();

	}

	delete [] Start;
	delete [] End;
	delete [] PASS1;
	delete [] PASS2;

	return 0;
}

int RunFQ2FQ_SEOUT ( Para_A24 * P2In,  vector<std::string>  & FilePath, const kc_c4x_t * h )
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAQQQ ;
	vector <string>  AAAIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE  > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAQQQ.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);


	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];


	bool  *PASS1 =new bool [BATCH_SIZE];

	igzstream INA (FilePath[0].c_str(),ifstream::in); // ifstream  + gz
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

	string ID_1 ,seq_1,temp_1,Quly_1;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+".fq";

	if (P2In->OUTGZ)
	{
		OUTIOGZ OUTHanDle ;
		outputFileAPE=OUT+".fq.gz";

		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());		
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.SEAA=0;

		int CountFQ=0;

		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j])
					{
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.SEAA++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				{
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					OUTHanDle.SEAA++;
				}
			}
			CountFQ=0;
		}

		cout<<"INFO: output SE read number is "<<OUTHanDle.SEAA<<"\n";

		OUTHanDle.OUTAPE.close();

	}
	else
	{
		OUTIO OUTHanDle ;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());		
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.SEAA=0;
		int CountFQ=0;

		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);
			getline(INA,temp_1);
			getline(INA,Quly_1);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAQQQ[CountFQ]=Quly_1;
			AAAIII[CountFQ]=ID_1;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}
				threads.clear();
				RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j])
					{
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
						OUTHanDle.SEAA++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j])
				{
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n+\n"<<AAAQQQ[j]<<"\n";
					OUTHanDle.SEAA++;
				}
			}
			CountFQ=0;
		}

		cout<<"INFO: output SE read number is "<<OUTHanDle.SEAA<<"\n";

		OUTHanDle.OUTAPE.close();

	}

	delete [] Start;
	delete [] End;
	delete [] PASS1;

	return 0;
}

int RunFA2FA_PEOUT ( Para_A24 * P2In,  vector<std::string>  & FilePath, const kc_c4x_t *h )
{


	std::ios::sync_with_stdio(false);
	std::cin.tie(0);


	vector <string>  AAASSS ;
	vector <string>  AAAIII ;


	vector <string>  BBBSSS ;
	vector <string>  BBBIII ;


	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if (BATCH_SIZE > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	BBBSSS.resize(BATCH_SIZE+2);
	BBBIII.resize(BATCH_SIZE+2);


	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];
	bool  *PASS1 =new bool [BATCH_SIZE+2];
	bool  *PASS2 =new bool [BATCH_SIZE+2];




	igzstream INA (FilePath[0].c_str(),ifstream::in); // ifstream  + gz
	igzstream INB (FilePath[1].c_str(),ifstream::in); // ifstream  + gz
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
	INB.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);


	string ID_1 ,seq_1;
	string ID_2 ,seq_2;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+"_pe_1.fa";
	string outputFileBPE=OUT+"_pe_2.fa";
	string outputFileASE=OUT+"_se_1.fa";
	string outputFileBSE=OUT+"_se_2.fa";

	if (P2In->OUTGZ)
	{
		OUTIOGZ OUTHanDle ;
		outputFileAPE=OUT+"_pe_1.fa.gz";
		outputFileBPE=OUT+"_pe_2.fa.gz";
		outputFileASE=OUT+"_se_1.fa.gz";
		outputFileBSE=OUT+"_se_2.fa.gz";

		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		int CountFQ=0;

		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);

			getline(INB,ID_2);
			getline(INB,seq_2);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAIII[CountFQ]=ID_1;

			BBBSSS[CountFQ]=seq_2;
			BBBIII[CountFQ]=ID_2;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j]  & PASS2[j] )
				{						
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.PE++;
				}
				else if (PASS1[j])
				{
					OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.SEAA++;
				}
				else if (PASS2[j])
				{
					OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.SEBB++;
				}
			}
			CountFQ=0;
		}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();

	}
	else
	{
		OUTIO OUTHanDle ;

		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTASE.open(outputFileASE.c_str());
		OUTHanDle.OUTBPE.open(outputFileBPE.c_str());
		OUTHanDle.OUTBSE.open(outputFileBSE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTASE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.OUTBSE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		OUTHanDle.PE=0;OUTHanDle.SEAA=0;OUTHanDle.PE=0;OUTHanDle.SEBB=0;
		int CountFQ=0;
		//		cerr<<"Start\t"<<CountFQ<<"\t"<<BATCH_SIZE<<"\t"<<n_thread<<endl;
		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);

			getline(INB,ID_2);
			getline(INB,seq_2);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAIII[CountFQ]=ID_1;

			BBBSSS[CountFQ]=seq_2;
			BBBIII[CountFQ]=ID_2;

			CountFQ++;
			//			cerr<<CountFQ<<"Start\t"<<CountFQ<<"\t"<<BATCH_SIZE<<"\t"<<n_thread<<endl;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j]  & PASS2[j] )
					{						
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.PE++;
					}
					else if (PASS1[j])
					{
						OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
					else if (PASS2[j])
					{
						OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
						OUTHanDle.SEBB++;
					}
				}
				CountFQ=0;
			}
		}
		
		CountFQ--;
		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				
				threads.push_back(std::thread(FilterLowHitPE,P2In ,std::ref(h),PASS1,PASS2,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS),std::ref(BBBSSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRPE(P2In,PASS1,PASS2,CountFQ,AAASSS,BBBSSS);
			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j]  & PASS2[j] )
				{						
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.OUTBPE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.PE++;
				}
				else if (PASS1[j])
				{
					OUTHanDle.OUTASE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.SEAA++;
				}
				else if (PASS2[j])
				{
					OUTHanDle.OUTBSE<< BBBIII[j]<<"\n"<<BBBSSS[j]<<"\n";
					OUTHanDle.SEBB++;
				}
			}
			CountFQ=0;
		}

		cout<<"INFO: output PE1 read number is "<<OUTHanDle.PE+OUTHanDle.SEAA<<"\n";
		cout<<"INFO: output PE2 read number is "<<OUTHanDle.PE+OUTHanDle.SEBB<<"\n";
		cout<<"INFO: paired PE1 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: paired PE2 read number is "<<OUTHanDle.PE<<"\n";
		cout<<"INFO: un-paired PE1 read number is "<<OUTHanDle.SEAA<<"\n";
		cout<<"INFO: un-paired PE2 read number is "<<OUTHanDle.SEBB<<"\n";

		OUTHanDle.OUTAPE.close();
		OUTHanDle.OUTASE.close();
		OUTHanDle.OUTBPE.close();
		OUTHanDle.OUTBSE.close();

	}

	delete [] Start;
	delete [] End;
	delete [] PASS1;
	delete [] PASS2;

	return 0;
}

int RunFA2FA_SEOUT ( Para_A24 * P2In,  vector<std::string>  & FilePath, const kc_c4x_t *h )
{

	std::ios::sync_with_stdio(false);
	std::cin.tie(0);

	vector <string>  AAASSS ;
	vector <string>  AAAIII ;

	int BinWind=VECMAX; 
	int BATCH_SIZE;
	BATCH_SIZE=BinWind*n_thread;

	if  (BATCH_SIZE  > (P2In->ReadNumber)) 
	{ 
		BinWind =(P2In->ReadNumber)/n_thread; 
		if (BinWind<2) {BinWind=2;} ;
		BATCH_SIZE=BinWind*n_thread;
	}

	AAASSS.resize(BATCH_SIZE+2);
	AAAIII.resize(BATCH_SIZE+2);

	std::vector<std::thread> threads;

	int * Start =new int [n_thread];
	int * End =new int [n_thread];

	bool  *PASS1 =new bool [BATCH_SIZE];

	igzstream INA (FilePath[0].c_str(),ifstream::in); // ifstream  + gz
	INA.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

	string ID_1 ,seq_1;
	string ID_2 ,seq_2;

	string OUT=(P2In->OutFq1);

	string outputFileAPE=OUT+".fa";

	if (P2In->OUTGZ)
	{
		OUTIOGZ OUTHanDle ;
		outputFileAPE=OUT+".fa.gz";
		OUTHanDle.SEAA=0;
		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);

		int CountFQ=0;

		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);

			if (ID_1.empty()) {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAIII[CountFQ]=ID_1;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();

				RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j])
					{
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}
			threads.clear();
			RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j])
				{
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.SEAA++;
				}
			}
			CountFQ=0;
		}

		cout<<"INFO: output SE read number is "<<OUTHanDle.SEAA<<"\n";

		OUTHanDle.OUTAPE.close();

	}
	else
	{
		OUTIO OUTHanDle ;

		OUTHanDle.OUTAPE.open(outputFileAPE.c_str());
		OUTHanDle.OUTAPE.rdbuf()->pubsetbuf(nullptr, BATCH_SIZE*1024);
		OUTHanDle.SEAA=0;
		int CountFQ=0;

		while(!INA.eof())
		{
			getline(INA,ID_1);
			getline(INA,seq_1);

			if (ID_1.empty())   {continue ; }
			AAASSS[CountFQ]=seq_1;
			AAAIII[CountFQ]=ID_1;

			CountFQ++;
			if (CountFQ==BATCH_SIZE)
			{
				for (int i = 0; i < n_thread; i++)
				{
					Start[i]=i*BinWind;
					End[i]=Start[i]+BinWind;
					threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
				}

				for (auto& thread : threads)
				{
					thread.join();
				}

				threads.clear();
				RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

				for (int j = 0; j < CountFQ; j++)
				{
					if (PASS1[j])
					{
						OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
						OUTHanDle.SEAA++;
					}
				}
				CountFQ=0;
			}
		}

		if (CountFQ!=0)
		{
			for (int i = 0; i < n_thread; i++)
			{
				Start[i]=i*BinWind;
				End[i]=Start[i]+BinWind;
				if (End[i]>CountFQ)  {  End[i]=CountFQ;  }
				threads.push_back(std::thread(FilterLowHitSE,P2In ,std::ref(h),PASS1,std::ref(Start[i]),std::ref(End[i]),std::ref(AAASSS)));
			}

			for (auto& thread : threads)
			{
				thread.join();
			}

			threads.clear();
			RmPCRSE(P2In,PASS1,CountFQ,AAASSS);

			for (int j = 0; j < CountFQ; j++)
			{
				if (PASS1[j])
				{
					OUTHanDle.OUTAPE<< AAAIII[j]<<"\n"<<AAASSS[j]<<"\n";
					OUTHanDle.SEAA++;
				}

			}
			CountFQ=0;
		}

		cout<<"INFO: output SE read number is "<<OUTHanDle.SEAA<<"\n";

		OUTHanDle.OUTAPE.close();

	}

	delete [] Start;
	delete [] End;
	delete [] PASS1;

	return 0;
}

//////////////////main///////////////////
int main (int argc, char *argv[ ])
{

	Para_A24 * P2In = new Para_A24;
	int InPESE=1;
	InPESE=parse_cmd_FqSplit(argc, argv, P2In );
	if(InPESE==0)
	{
		delete P2In ;
		return 0 ;
	}

	string path=(P2In->OutFq1);
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (ext == "gz")
	{
		(P2In->OutFq1)=path.substr(0,path.rfind('.') ==string::npos ? path.length() : path.rfind('.'));
	}

	path=(P2In->OutFq1);
	ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	if (ext == "fa"  ||   ext == "fq")
	{
		(P2In->OutFq1)=path.substr(0,path.rfind('.') ==string::npos ? path.length() : path.rfind('.'));
	}

	if  (InPESE==2)
	{
		(P2In->InFq1)=P2In->InSeFq;
	}

	P2In->LowQint=GetShiftQ((P2In->InFq1),P2In); // (Phred) 33 or 64
	P2In->MinBaseQ=(P2In->MinBaseQ)+(P2In->LowQint);
	P2In->AverQ=(P2In->AverQ)+(P2In->LowQint);

	if ((P2In->HalfReadLength) < ((P2In->Kmer)+1))
	{
		P2In->HalfReadLength=(P2In->ReadLength)/2;
		if ((P2In->HalfReadLength)<((P2In->Kmer)+1) )
		{
			(P2In->HalfReadLength)=(P2In->Kmer)+1;
		}

		if  (((P2In->MinBaseQ)<2) && ((P2In->AverQ)<2)  && ((P2In->N_Number)==0))
		{
			P2In->FILTER=false ;
		}
	}

	for (int i=0; i<256; i++)
	{
		NArry[i]=false;
		LowArry[i]=true;
	}
	NArry['N']=true ; NArry['n']=true ;

	for (int i=0; i<(P2In->MinBaseQ); i++)
	{
		LowArry[i]=false;
	}

	if ((P2In->MinCount)==1)
	{
		RunM1_Work(P2In, InPESE);
		delete P2In ;
		return (0);
	}

	vector<std::string> FilePath;

	if ((InPESE==1) && ((P2In->LowQint)!=0)) // PE FQ
	{
		P2In->ReadNumber=(P2In->ReadNumber)/2;
		RunFQFilterPE(P2In,FilePath);
	}
	else if ((InPESE==2) && ((P2In->LowQint)!=0))  // SE FQ
	{
		RunFQFilterSE(P2In,FilePath);
	}
	else if ((InPESE==1) && ((P2In->LowQint)==0)) // PE FA
	{
		P2In->ReadNumber=(P2In->ReadNumber)/2;
		RunFAFilterPE(P2In,FilePath);
	}
	else if ((InPESE==2) && ((P2In->LowQint)==0)) // SE FA
	{
		RunFAFilterSE(P2In,FilePath);
	}

	uint64_t hash_size = 100000000; // Initial size of hash.  100M 
	if  (n_thread>10)
	{
		hash_size=int(n_thread*hash_size/10);
	}

	kc_c4x_t *h;
	int p=KC_BITS ;  // Minimum length of counting field

	// create the hash
	h = count_file(FilePath, P2In->Kmer, P2In->Windows, p, hash_size, n_thread);

	GetMinCount(P2In,h);

	if ((InPESE==1) && ((P2In->LowQint)!=0)) // PE FQ
	{
		if (P2In->OutFa)
		{
			RunFA2FA_PEOUT(P2In, FilePath, h);
		}
		else
		{
			RunFQ2FQ_PEOUT(P2In, FilePath, h);
		}
	}
	else if  ((InPESE==2) && ((P2In->LowQint)!=0))  // SE FQ
	{
		if (P2In->OutFa)
		{
			RunFA2FA_SEOUT(P2In, FilePath, h);
		}
		else
		{
			RunFQ2FQ_SEOUT(P2In, FilePath, h);
		}
	}
	else if ((InPESE==1) && ((P2In->LowQint)==0)) // PE FA
	{
		RunFA2FA_PEOUT(P2In, FilePath, h);
	}
	else if  ((InPESE==2) && ((P2In->LowQint)==0)) // SE FA
	{
		RunFA2FA_SEOUT(P2In, FilePath, h);
	}

	for (auto& file : FilePath)
	{
		remove(file.c_str());
	}

	for (int i = 0; i < 1<<p; ++i)
	{
		kc_c4_destroy(h->h[i]);
	}
	free(h->h); free(h);

	delete P2In ;
	return (0);
}

#endif  //
