/* The MIT License

   Copyright (c) 2023- by 
   Huiyang Yu,
   Weiming He,
   Chunmei Shi.

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

#ifndef bamCov_H_
#define bamCov_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <list>
#include <map>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include <unordered_map>
#include <thread>
#include <algorithm>
#include "./comm.h"
#include "./DataClass.h"
#include "./zlib.h"
#include "./sam.h"
#include "./hts.h"
#include "./cram.h"
#include "./kseq.h"
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

using namespace std;
typedef unsigned long long ubit64_t;

void bamCov_help()
{
	cout<<""
		"Usage: pandepth -i in.bam [-g gene.gff | -b region.bed] -o outPrefix\n"
		" Input/Output options:\n"
		"   -i    <str>     input of sam/bam/cram/paf or #.list file\n"
		"   -o    <str>     prefix of output file\n"
		" Target options:\n"
		"   -g    <str>     input gff/gtf file for gene region\n"
		"   -f    <str>     gff/gtf feature type to parse, CDS or exon [CDS]\n"
		"   -b    <str>     input bed file for list of regions\n"
		"   -w    <int>     windows size (bp)\n"
		"   -a              output all the site depth\n"
		" Filter options:\n"
		"   -q    <int>     min mapping quality [0]\n"
		"   -d    <int>     min site depth for statistics [1]\n"
		"   -x    <int>     exclude reads with any of the bits in FLAG set [1796]\n"
		" Other options:\n"
		"   -t    <int>     number of threads [3]\n"
		"   -r    <str>     reference genome file for cram decode or GC parse\n"
		"   -c              enable the calculation of GC content (requires -r)\n"
		"   -h              show this help [v2.26]\n"
		"\n";
}

int bamCov_help01(int argc, char **argv , In3str1v * paraFA04   )
{
	if (argc <=1 ) {bamCov_help();return 0;}
	int file_count=0;
	int Infile_count=0;
	for(int i = 1; i < argc ; i++)
	{
		if(argv[i][0] != '-')
		{
			cerr << "Error: Command option error! Please check the provided options." << endl;
			return 0;
		}
		string flag=argv[i] ;
		flag=replace_all(flag,"-","");

		if (flag  == "i" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			string ext =A.substr(A.rfind('.') ==string::npos ? A.length() : A.rfind('.') + 1);
			if ((ext == "list")  ||  (ext == "List" ))
			{
				ReadList(A, Infile_count ,(paraFA04->InStr1List) ) ;
			}
			else
			{
				Infile_count++;
				(paraFA04->InStr1List).push_back(A);
				paraFA04->InStr1=argv[i];
			}
		}
		else if (flag  ==  "o" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr3=argv[i];
		}
		else if (flag  ==  "c" )
		{
			paraFA04->gc=true;
		}
		else if (flag  ==  "a" )
		{
			paraFA04->SiteOutPut=true;
		}
		else if (flag  ==  "r" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->reference=argv[i];
		}
		else if (flag  ==  "f" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->CDS=argv[i];
		}
		else if (flag  ==  "x" )
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->flags =atoi(argv[i]);
		}

		else if (flag  =="g")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InStr2=argv[i];
			igzstream INGFF (argv[i],ifstream::in);
			if (INGFF.fail())
			{
				cerr << "Error: Failed to open the file: "<<argv[i]<<endl;
				return  0;
			}

			string tmp ;
			for (int A=1 ; A<168 && (!INGFF.eof())  ; A++ )
			{
				getline(INGFF,tmp);			
				if (tmp.length()<2)
				{
					continue ;
				}
				if (tmp[0] == '#')
				{
					continue ;
				}
				if (tmp.find("Parent")!=string::npos)
				{
					paraFA04->InInt2=1 ;
				}
				else if (tmp.find("transcript_id")!=string::npos)
				{
					paraFA04->InInt2=2 ;
				}
			}

			INGFF.close();

			if ((paraFA04->InInt2)==0)
			{
				cerr<<"Error: The format of the input GFF/GTF file is incorrect. Please check the file format: "<<argv[i]<<endl;
				return 0;
			}
		}
		else if (flag == "b")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			string A=argv[i];
			file_count++;
			(paraFA04->List).push_back(A);
			(paraFA04->InInt2)=3;
		}
		else if (flag == "t")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->CPU=atoi(argv[i]);
		}
		else if (flag == "w")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->WinSize=atoi(argv[i]);
			if ( ( paraFA04->WinSize) <1)  
			{
				cerr<<"Warning: -w should >= 1, set to 1\n";
				paraFA04->WinSize=1;
			}
		}
		else if (flag == "q")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->InInt=atoi(argv[i]);
		}
		else if (flag == "s")
		{
			paraFA04->TF=false ;
		}

		else if (flag == "d")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->minDep=atoi(argv[i]);
			if (paraFA04->minDep < 1) paraFA04->minDep = 1;
		}

		else if (flag == "help"  || flag == "h")
		{
			bamCov_help();return 0;
		}
		else
		{
			cerr << "Error UnKnow argument -"<<flag<<endl;
			return 0;
		}
	}
	//cerr<<Infile_count<<endl;
	if (Infile_count>0)
	{
		paraFA04->InStr1=(paraFA04->InStr1List)[0];
	}
	if (  (paraFA04->InStr1).empty()  || (paraFA04->InStr3).empty() )
	{
		cerr<< "Error: lack argument -i or -o "<<endl ;
		return 0;
	}

	if (  (paraFA04->InStr2).empty()  && file_count==0 )
	{
		//	cerr<< "lack argument for the must"<<endl ;
		//	return 0;
	}
	else if  ( file_count!=0  &&  (paraFA04->InStr2).empty() )
	{
		string ListName=(paraFA04->List)[0];
		paraFA04->InStr2=ListName;
		(paraFA04->List).clear();

		igzstream LISTTT (ListName.c_str(),ifstream::in); // igzstream
		if (!LISTTT.good())
		{
			cerr << "Error: Failed to open BED file: "<<ListName<<endl;
			return  0;
		}

		string  tmp,line;
		getline(LISTTT,tmp);
		getline(LISTTT,line);
		LISTTT.close();

		vector<string> Temp1;
		vector<string> Temp2;
		split(tmp,Temp1," \t");

		split(line,Temp2," \t");
		if 	(Temp1.size()==4  ||  Temp2.size()==4)
		{
			paraFA04->InInt2=4;
		}
	}

	(paraFA04->InStr3)=add_Asuffix(paraFA04->InStr3);
	return Infile_count ;
}

void  StatChrDepthWin ( unsigned int *depth,  map <string,GeneInfo> &  GeneStat,int MeMStart, int MeMEnd ,int ShiftPosition, int minDep)
{
	for (auto iter = GeneStat.begin(); iter != GeneStat.end(); ++iter)
	{
		if ( (MeMEnd-MeMStart)!=0)
		{
			if ( (iter->second).GeneStart >=MeMEnd ) { continue;}
			if ( (iter->second).GeneEnd <MeMStart ) { continue;}
		}
		else
		{
			if ( (iter->second).GeneStart >MeMEnd ) { continue;}
			if ( (iter->second).GeneEnd <MeMStart ) { continue;}
		}
		int Size=(iter->second).CDSList.size();
		//cerr<<"\tww\t"<<iter->first<<"\t"<< (iter->second).GeneStart<<"\t"<< (iter->second).GeneEnd<<endl;
		for (int tt=0  ; tt< Size ; tt++)
		{
			int ii=((iter->second).CDSList[tt].first-1)-ShiftPosition;
			int End=((iter->second).CDSList[tt].second)-ShiftPosition;


			for (  ; ii<End ; ii++)
			{
				if ( depth[ii] >= minDep)
				{
					(iter->second).GeneCover++;
					(iter->second).GeneDepth+=depth[ii];
				}
			}
		}
	}
}

void  StatChrDepthLowMEM (SiteInfo  * depth  ,  map <string,GeneInfo>  &  GeneStat, int minDep)
{
	for (auto iter = GeneStat.begin(); iter != GeneStat.end(); ++iter)
	{
		int Size=(iter->second).CDSList.size();
		for (int tt=0  ; tt< Size ; tt++)
		{
			int ii=(iter->second).CDSList[tt].first-1;
			int End=(iter->second).CDSList[tt].second ;
			for (  ; ii<End ; ii++)
			{
				if ( (depth[ii].Depth) >= minDep)
				{
					(iter->second).GeneCover++;
					(iter->second).GeneDepth+=(depth[ii].Depth);
				}
			}
		}
	}
}

void ProDealChrBambaiOUTSite ( string  & BamPath , In3str1v * paraFA04 , SiteInfo **  depth , map <int,map <int,int> >  & RegionMerger ,  vector <int> &  ChrNumVer , map <int,map <string,GeneInfo> > & GeneData , int &  numThreads )
{

	htsFile *fphts;
	sam_hdr_t *headerAA;
	int64_t rcnt;

	int n=1;
	int i=0;
	hts_idx_t *idx=NULL;
	int64_t *cnt ;
	int ret ;

	int min_mapQ= (paraFA04->InInt);
	fphts = hts_open(BamPath.c_str(), "r");

	hts_set_log_level(HTS_LOG_OFF);
	//
	if(fphts->format.format == htsExactFormat::cram)
	{
		const char* ref_file = (paraFA04->reference).c_str();
		hts_set_fai_filename(fphts, ref_file);
		hts_set_opt(fphts, CRAM_OPT_DECODE_MD, 0);
		hts_set_opt(fphts, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	}

	if (fphts)
	{
		idx=sam_index_load(fphts,BamPath.c_str());
	}

	if (fphts == 0 || idx == 0 ) 
	{
		cerr<<"Error: Failed to open the index file or BAM/CRAM file: "<<BamPath<<endl; 
		return ;
	}
	//
	hts_set_opt(fphts, HTS_OPT_NTHREADS, numThreads);
	hts_set_opt(fphts, CRAM_OPT_DECODE_MD, 0);
	hts_set_opt(fphts, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	headerAA = sam_hdr_read(fphts);

	if (headerAA == NULL) { cerr<<"Error: Failed to read the header for the BAM/CRAM file: "<<BamPath<<endl; return  ;}

	bam1_t *aln = bam_init1();
	uint32_t flags = (paraFA04->flags);
	hts_itr_t *iter;
	uint32_t *cigar;

	int AAA=ChrNumVer.size();
	map <int,map <int,int> > :: iterator  RegionIt ;
	map <int,int> :: iterator MapSSEE;

	for (int po=0; po<AAA ; po++)
	{
		int ChrNum=ChrNumVer[po];
		int ChrLen=(headerAA->target_len[ChrNum]);
		string ChrName=headerAA->target_name[ChrNum];

		RegionIt =  RegionMerger.find(ChrNum);
		MapSSEE =(RegionIt->second).begin();
		int NumberRegionSize=(RegionIt->second).size();

		char* CharMap = new char[NumberRegionSize*128];
		char ** RegionArry = new  char * [NumberRegionSize];

		uint64_t CountRegion = 0;


		for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++)
		{
			long beg=(MapSSEE->first)-1;
			if  (beg<1) {beg=1;}
			long end=(MapSSEE->second)+1;
			if (end>ChrLen) {end=ChrLen;}
			RegionArry[CountRegion++] = CharMap;
			//cerr<<ChrName<<"\t"<<beg<<"\t"<<endl;
			CharMap += (sprintf(CharMap, "%s:%lu-%lu", ChrName.c_str(), beg, end)+1);
		}
		rcnt = CountRegion;
		iter = sam_itr_regarray(idx, headerAA, RegionArry, rcnt);

		int32_t endTmp;
		int32_t StartRead;
		while ((ret = sam_itr_next(fphts, iter, aln)) >= 0)
		{
			if ( aln->core.flag & flags )  { continue; }
			if ( (aln->core).qual < (paraFA04->InInt) )	{	continue ;}
			cigar = bam_get_cigar(aln);
			StartRead=((aln->core).pos);
			for(int i=0; i < aln->core.n_cigar;++i)
			{
				int cig=bam_cigar_op(cigar[i]);
				int ncig = bam_cigar_oplen(cigar[i]);
				switch (cig)
				{
					case 0:
					case 7:
					case 8:
						endTmp=StartRead+ncig;
						for (  ; StartRead<endTmp;StartRead++)
						{
							(depth[ChrNum][StartRead]).Depth++;
						}
						break;
					case 2:
					case 3:
						StartRead=StartRead+ncig;
						break;
				}
			}
		}


		auto GeneDataIT=GeneData.find(ChrNum);


		StatChrDepthLowMEM( depth[ChrNum]  , GeneDataIT->second, paraFA04->minDep);
		delete [] RegionArry ;

	}


	//delete CharMap ;


	bam_destroy1(aln);
	if (idx) hts_idx_destroy(idx);
	if (iter) sam_itr_destroy(iter);
	if (headerAA) sam_hdr_destroy(headerAA);




}


void ProDealChrBambaiALLSite ( string  & BamPath , In3str1v * paraFA04 , SiteInfo **  depth , map <int,map <int,int> >  & RegionMerger ,  vector <int> &  ChrNumVer , map <int,map <string,GeneInfo> > & GeneData , int &  numThreads )
{

	htsFile *fphts;
	sam_hdr_t *headerAA;
	int64_t rcnt;

	int n=1;
	int i=0;
	hts_idx_t *idx=NULL;
	int64_t *cnt ;
	int ret ;

	int min_mapQ= (paraFA04->InInt);
	fphts = hts_open(BamPath.c_str(), "r");

	hts_set_log_level(HTS_LOG_OFF);
	//
	if(fphts->format.format == htsExactFormat::cram)
	{
		const char* ref_file = (paraFA04->reference).c_str();
		hts_set_fai_filename(fphts, ref_file);
		hts_set_opt(fphts, CRAM_OPT_DECODE_MD, 0);
		hts_set_opt(fphts, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	}

	if (fphts)
	{
		idx=sam_index_load(fphts,BamPath.c_str());
	}

	if (fphts == 0 || idx == 0 ) 
	{
		cerr<<"Error: Failed to open the index file or BAM/CRAM file: "<<BamPath<<endl; 
		return ;
	}
	//
	hts_set_opt(fphts, HTS_OPT_NTHREADS, numThreads);
	hts_set_opt(fphts, CRAM_OPT_DECODE_MD, 0);
	hts_set_opt(fphts, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	headerAA = sam_hdr_read(fphts);

	if (headerAA == NULL) { cerr<<"Error: Failed to read the header for the BAM/CRAM file: "<<BamPath<<endl; return  ;}

	bam1_t *aln = bam_init1();
	uint32_t flags = (paraFA04->flags);
	hts_itr_t *iter;
	uint32_t *cigar;

	int AAA=ChrNumVer.size();
	map <int,map <int,int> > :: iterator  RegionIt ;
	map <int,int> :: iterator MapSSEE;

	for (int po=0; po<AAA ; po++)
	{
		int ChrNum=ChrNumVer[po];
		int ChrLen=(headerAA->target_len[ChrNum]);
		string ChrName=headerAA->target_name[ChrNum];

		RegionIt =  RegionMerger.find(ChrNum);
		MapSSEE =(RegionIt->second).begin();
		int NumberRegionSize=(RegionIt->second).size();

		char* CharMap = new char[NumberRegionSize*128];
		char ** RegionArry = new  char * [NumberRegionSize];

		uint64_t CountRegion = 0;


		for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++)
		{
			long beg=(MapSSEE->first)-1;
			if  (beg<1) {beg=1;}
			long end=(MapSSEE->second)+1;
			if (end>ChrLen) {end=ChrLen;}
			RegionArry[CountRegion++] = CharMap;
			CharMap += (sprintf(CharMap, "%s:%lu-%lu", ChrName.c_str(), beg, end)+1);
		}
		rcnt = CountRegion;
		iter = sam_itr_regarray(idx, headerAA, RegionArry, rcnt);

		int32_t endTmp;
		int32_t StartRead;
		while ((ret = sam_itr_next(fphts, iter, aln)) >= 0)
		{
			if ( aln->core.flag & flags )  { continue; }
			if ( (aln->core).qual < (paraFA04->InInt) )	{	continue ;}
			cigar = bam_get_cigar(aln);
			StartRead=((aln->core).pos);
			for(int i=0; i < aln->core.n_cigar;++i)
			{
				int cig=bam_cigar_op(cigar[i]);
				int ncig = bam_cigar_oplen(cigar[i]);
				switch (cig)
				{
					case 0:
					case 7:
					case 8:
						endTmp=StartRead+ncig;
						for (  ; StartRead<endTmp;StartRead++)
						{
							(depth[ChrNum][StartRead]).Depth++;
						}
						break;
					case 2:
					case 3:
						StartRead=StartRead+ncig;
						break;
				}
			}
		}


		delete [] RegionArry ;

	}




	bam_destroy1(aln);
	if (idx) hts_idx_destroy(idx);
	if (iter) sam_itr_destroy(iter);
	if (headerAA) sam_hdr_destroy(headerAA);




}


void ProDealChrBambai ( string  & BamPath , In3str1v * paraFA04   ,  map <int,map <int,int> >  & RegionMerger , vector <int> &  ChrNumVer ,  map <int,map <string,GeneInfo> > & GeneData , int &  numThreads  )
{

	htsFile *fphts;
	sam_hdr_t *headerAA;
	int64_t rcnt;


	int n=1;
	int i=0;
	hts_idx_t *idx=NULL;
	int64_t *cnt ;
	int ret ;

	int min_mapQ= (paraFA04->InInt);
	fphts = hts_open(BamPath.c_str(), "r");

	hts_set_log_level(HTS_LOG_OFF);
	//
	if(fphts->format.format == htsExactFormat::cram)
	{
		const char* ref_file = (paraFA04->reference).c_str();
		hts_set_fai_filename(fphts, ref_file);
		hts_set_opt(fphts, CRAM_OPT_DECODE_MD, 0);
		hts_set_opt(fphts, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	}
	//

	if (fphts)
	{
		//  idx=sam_index_load3(fphts, BamPath.c_str(), NULL, HTS_IDX_SILENT_FAIL);
		idx=sam_index_load(fphts,BamPath.c_str());
	}

	if (fphts == 0 || idx == 0 ) 
	{
		cerr<<"Error: Failed to open the index file or BAM/CRAM file: "<<BamPath<<endl; 
		return ;
	}
	//
	hts_set_opt(fphts, HTS_OPT_NTHREADS, numThreads);
	hts_set_opt(fphts, CRAM_OPT_DECODE_MD, 0);
	hts_set_opt(fphts, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	headerAA = sam_hdr_read(fphts);

	if (headerAA == NULL) { cerr<<"Error: Failed to read the header for the BAM/CRAM file: "<<BamPath<<endl; return  ;}

	bam1_t *aln = bam_init1();
	uint32_t flags = (paraFA04->flags);
	hts_itr_t *iter;
	uint32_t *cigar;

	int AAA=ChrNumVer.size();
	map <int,map <int,int> > :: iterator  RegionIt ;
	map <int,int> :: iterator MapSSEE ;
	map <int,int> :: iterator MapSSEEV2 ;

	//int  MeMBinWindows=int (ChrLen/40) ; 
	int  MeMBinWindows=10000000; 
	//if  (MeMBinWindows<20000000) {MeMBinWindows=20000000;}
	//int  MeMBinWindowsEdge=int(MeMBinWindows/10);
	int  MeMBinWindowsEdge=10000000;

	for (int po=0; po<AAA ; po++)
	{
		int ChrNum=ChrNumVer[po];
		int ChrLen=(headerAA->target_len[ChrNum]);
		string ChrName=headerAA->target_name[ChrNum];

		RegionIt =  RegionMerger.find(ChrNum) ;
		MapSSEE =(RegionIt->second).begin() ;
		int NumberRegionSize=(RegionIt->second).size();

		char* CharMap = new char[NumberRegionSize*128];
		char ** RegionArry = new  char * [NumberRegionSize];

		uint64_t CountRegion = 0;
		int MeMStart=(MapSSEE->first); if (MeMStart<1) {MeMStart=1;}
		int MeMEnd=MeMStart+MeMBinWindows-1; if (MeMEnd>ChrLen) {MeMEnd=ChrLen;}

		for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++)
		{
			long beg=(MapSSEE->first)-1;
			if  (beg<1) {beg=1;}
			long end=(MapSSEE->second)+1;
			if (end>ChrLen) {end=ChrLen;}
			RegionArry[CountRegion++] = CharMap;
			//cerr<<ChrName<<"\t"<<beg<<"\t"<<end<<"\n";
			CharMap += (sprintf(CharMap, "%s:%lu-%lu", ChrName.c_str(), beg, end)+1);

			MapSSEEV2=MapSSEE;
			MapSSEEV2++;

			if ( end>=MeMEnd  ||   (MapSSEEV2==(RegionIt->second).end()) )
			{
				MeMEnd=end;

				int ThisWinDowLen=MeMEnd-MeMStart+1;
				ThisWinDowLen+=(2*MeMBinWindowsEdge);
				unsigned int *depth = new unsigned int [ThisWinDowLen];
				for (int  ccv=0 ; ccv<ThisWinDowLen  ; ccv++)
				{
					depth[ccv]=0;
				}

				int  ShiftPosition=MeMStart-MeMBinWindowsEdge-1;

				rcnt = CountRegion;
				iter = sam_itr_regarray(idx, headerAA, RegionArry, rcnt);



				while ((ret = sam_itr_next(fphts, iter, aln)) >= 0)
				{
					if ( aln->core.flag & flags )  { continue; }
					if ( (aln->core).qual < (paraFA04->InInt) )	{	continue ;	}

					cigar = bam_get_cigar(aln);
					int32_t StartRead=((aln->core).pos)-ShiftPosition;
					int32_t endTmp;
			//		cerr<<(aln->core).pos<<"\the"<<endl;
					for(int i=0; i < aln->core.n_cigar;++i)
					{
						int cig=bam_cigar_op(cigar[i]);
						int ncig = bam_cigar_oplen(cigar[i]);
						switch (cig)
						{
							case 0:
							case 7:
							case 8:
								endTmp=StartRead+ncig;
								for (  ; StartRead<endTmp;StartRead++)
								{
									depth[StartRead]++;
								}
								break;
							case 2:
							case 3:
								StartRead=StartRead+ncig;
								break;
						}

					}
				}


				auto GeneDataIT=GeneData.find(ChrNum);
				//cerr<<ChrNum<<"\t"<<MeMStart<<"\t"<<MeMEnd<<endl;
				StatChrDepthWin( depth  , GeneDataIT->second ,MeMStart,MeMEnd,ShiftPosition, paraFA04->minDep);

				delete [] depth;




				if (MapSSEEV2!=(RegionIt->second).end())
				{					
					MeMStart=MapSSEEV2->first;
					if (MeMStart-150 > MeMEnd )
					{
						MeMStart=MeMStart-150;
					}
				}
				MeMEnd=MeMStart+MeMBinWindows; if (MeMEnd>ChrLen) {MeMEnd=ChrLen;}
				CountRegion=0;
				CharMap=RegionArry[0];
			}

		}





		delete [] RegionArry ;
		//		delete CharMap ;

	}

	bam_destroy1(aln);
	if (idx) hts_idx_destroy(idx);
	if (iter) sam_itr_destroy(iter);
	if (headerAA) sam_hdr_destroy(headerAA);
}




bool parcigar(const std::string& IncG, std::vector<char>& charVec, std::vector<int>& intVec)
{
	size_t pos = IncG.find("cg:Z:");
	if (pos == std::string::npos) {     return false;    }

	std::string cgValue = IncG.substr(pos + 5);  // ��ȡ cg:Z: ����Ĳ���

	size_t length = cgValue.length();
	size_t start = 0, end = 0;

	while (start < length && end < length) 
	{
		while (end < length && isdigit(cgValue[end]))
		{
			end++;
		}

		int num = std::stoi(cgValue.substr(start, end - start));
		char operation = cgValue[end];
		charVec.push_back(operation);
		intVec.push_back(num);
		end++;
		start = end;
	}

	return true;
}





int findFirstStringIndex(const std::vector<std::string>& vec) 
{
	for (int i = 0; i < vec.size(); i++)
	{
		if (vec[i].substr(0, 5) == "cg:Z:") 
		{
			return i;

		}
	}
	return -1;
}


int paf_main(In3str1v *paraFA04 )
{

	string  BamPath=(paraFA04->InStr1);
	if (BamPath.length()<=0)  
	{
		cerr<<"Error: Failed to open the PAF file: "<<BamPath<<endl;
		return 1;
	}

	map <int,string>  RefBase ;
	bool  RefIn=false;


	map <string,int> Chr2IntMap; 
	map <int,int> target_len; 
	map <int,string> target_name;

	int ID=0;
	string chr;
	int chrlength;
	if (!(paraFA04->reference).empty())
	{

		gzFile fpRef;
		kseq_t *seq;
		int l;
		fpRef = gzopen((paraFA04->reference).c_str(), "r");
		seq = kseq_init(fpRef);
		RefIn=true;

		if ((paraFA04->gc)==true)
		{

			while ((l = kseq_read(seq)) >= 0 )
			{
				string chr=(seq->name.s);
				Chr2IntMap[chr]=ID;
				target_name[ID]=chr;
				target_len[ID]=seq->seq.l;
				string seqBB=seq->seq.s;
				RefBase.insert( map <int,string>  :: value_type (ID,seqBB));
				ID++;
			}
		}
		else
		{

			while ((l = kseq_read(seq)) >= 0 )
			{
				string chr=(seq->name.s);
				Chr2IntMap[chr]=ID;
				target_name[ID]=chr;
				target_len[ID]=seq->seq.l;
				ID++;
			}
		}
		kseq_destroy(seq); gzclose(fpRef);
	}
	else
	{
		if ((paraFA04->gc)==true)
		{
			cerr<< "Error: lack reference sequence (-r) for GC parse"<<endl ;
			return 0;
		}

		igzstream IN (BamPath.c_str(),ifstream::in);
		if(!IN.good())
		{
			cerr << "open PAF file error: "<<BamPath<<endl;
			return 1;
		}
		string tmp1,tmp2,tmp3,tmp4,tmp5;
		while(!IN.eof())
		{
			string  line ,chr_id ;
			int position ;
			getline(IN,line);
			if (line.length()<=0)  { continue  ; }
			istringstream isone (line,istringstream::in);
			isone>>tmp1>>tmp2>>tmp3>>tmp4>>tmp5>>chr>>chrlength;
			if (Chr2IntMap.find(chr)==Chr2IntMap.end())
			{
				Chr2IntMap[chr]=ID;
				target_name[ID]=chr;
				target_len[ID]=chrlength;
				ID++;
			}
		}
		IN.close();
	}

	int n_targets=ID;

	int GCGCArry[256] = {0};
	GCGCArry['C']=1; GCGCArry['c']=1; GCGCArry['G']=1; GCGCArry['g']=1; 
	GCGCArry['a']=0; GCGCArry['A']=0; GCGCArry['t']=0; GCGCArry['T']=0;

	map <string,int> :: iterator  MapItChr2Int ;
	unordered_map <string,int> :: iterator  UnMapItChr2Int ;
	int Start ; int End ;  string ChrName ;
	string GeneID ;
	map <int,map <string,GeneInfo> > GeneData;
	map <int,map <string,GeneInfo> >  :: iterator  GeneDataIT;
	map <string,GeneInfo>   :: iterator  GeneInfoIT  ;

	if  ((paraFA04->InInt2)!=0)
	{
		igzstream LIST ((paraFA04->InStr2).c_str(),ifstream::in); // igzstream
		if (!LIST.good())
		{
			cerr << "Error: Can't open the GFF/GTF File: "<<(paraFA04->InStr1)<<endl;
			return  1;
		}

		if ((paraFA04->InInt2)==1)
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				istringstream isone (line,istringstream::in);
				string flag , CDS , ZhengFu ,geneID ;;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End >>flag>>ZhengFu>>flag>>geneID;;

				vector<string> inf;
				vector<string> Temp;
				split(geneID,inf,",;");
				split(inf[0],Temp,"=");
				string GeneID= Temp[Temp.size()-1];
				for ( int jk=1;  jk<inf.size() ; jk++)
				{
					vector<string> Temptmp2;
					split(inf[jk],Temptmp2,"=");
					if (Temptmp2[0] == "Parent")
					{
						GeneID= Temptmp2[Temptmp2.size()-1];
					}
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{
					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];
							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}

				}
			}

		}
		else if ((paraFA04->InInt2)==2)
		{
			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				line=replace_all(line,"\"","");
				line=replace_all(line,";","");

				istringstream isone (line,istringstream::in);
				string  flag , CDS , ZhengFu ,geneID;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End >>flag>>ZhengFu>>flag  ;


				vector<string> inf;
				split(line,inf,"\t ");
				string GeneID=inf[9];
				for ( int jk=8;  jk<inf.size() ; jk+=2)
				{

					if (inf[0] == "transcript_id")
					{
						GeneID= inf[jk+1];
					}
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{
					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}


				}
			}
		}
		else if ((paraFA04->InInt2)==3)
		{
			string StartStr , EndStr ;
			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);
				if (line.length()<=0)  {continue;}
				if (line[0] == '#')  { continue;}
				istringstream isone (line,istringstream::in);
				isone>>ChrName>>StartStr>>EndStr;
				GeneID=ChrName+"_"+StartStr+"_"+EndStr ;
				Start=atoi(StartStr.c_str()); End=atoi(EndStr.c_str());
				if  (Start > End )
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
					continue;
				}
				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{

					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}
				}
			}
		}  // end (paraFA04->InInt2)==3

		else if ((paraFA04->InInt2)==4)
		{
			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);
				if (line.length()<=0)  {continue;}
				if (line[0] == '#')  { continue;}
				istringstream isone (line,istringstream::in);
				isone>>ChrName>>Start>>End>>GeneID;
				if  (Start > End )
				{
					cerr<<line<<"Warning: This region may be incorrect. \n"<<endl;
					continue;
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{
					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn)
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}
				}
			}
		}
		else
		{

		}

		LIST.close();
	}

	map <int,map <int,int> > RegionMerger;
	map <int,map <int,int> > :: iterator MergerIt;
	map <int,int> :: iterator  MergerMapSSEE;

	for ( GeneDataIT = GeneData.begin(); GeneDataIT != GeneData.end(); ++GeneDataIT)
	{
		int ChrInt = GeneDataIT ->first ;
		map <int,int>  RegionAAA;
		map <int,int> :: iterator MapSSEE ;
		for (GeneInfoIT=(GeneDataIT -> second).begin() ; GeneInfoIT!=(GeneDataIT -> second).end()  ;  GeneInfoIT ++)
		{
			Start=(GeneInfoIT->second).GeneStart;
			End=(GeneInfoIT->second).GeneEnd;
			MapSSEE=RegionAAA.find(Start);
			if (MapSSEE==RegionAAA.end())
			{
				RegionAAA.insert(map <int,int>  :: value_type(Start,End)) ;
			}
			else
			{
				if (End >(MapSSEE->second))
				{
					MapSSEE->second=End;
				}
			}
		}

		MapSSEE=RegionAAA.begin() ;
		Start=MapSSEE->first;
		End=MapSSEE->second;

		map <int,int> Start2End ;
		Start2End[Start]=End;
		RegionMerger.insert(map <int,map <int,int> > ::value_type(ChrInt,Start2End));
		MergerIt=RegionMerger.find(ChrInt);
		MapSSEE++;

		for(  ; MapSSEE!=RegionAAA.end() ; MapSSEE++ )
		{
			if (  (MapSSEE->first) >  End    )
			{
				Start=MapSSEE->first;
				End=MapSSEE->second;
				(MergerIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
			}
			else if  (  (MapSSEE->second) >  End    )
			{
				MergerMapSSEE=(MergerIt->second).end(); MergerMapSSEE--;
				MergerMapSSEE->second= (MapSSEE->second) ;
				End=MapSSEE->second;
			}
		}
	}

	// check region
	if (RegionMerger.empty())
	{
		int  MeMBinWindows=10000000;
		if ((paraFA04->WinSize) ==0 )
		{
			//cout<<"Warning: GFF/GTF or BED was not provided, the total chromosome length will be used as the parsing region.\n";
			(paraFA04->InInt2)=0;
		}
		else if ((paraFA04->WinSize) <150 )
		{
			(paraFA04->InInt2)=6;
		}
		else
		{
			(paraFA04->InInt2)=5;
			MeMBinWindows=(paraFA04->WinSize);
		}

		for(int i = 0; i < (n_targets); i++)
		{
			Start=1;
			End=2;

			for (Start=1;  End <= (target_len[i]) ; Start+=MeMBinWindows )
			{
				End=Start+MeMBinWindows-1;
				if (End>(target_len[i]))
				{
					End=(target_len[i]);
				}

				string GeneID=target_name[i]+Int2Str(Start);
				map <string,GeneInfo> TmpGene;
				GeneInfo  GeneInfoTmp ;
				GeneInfoTmp.GeneStart=Start;
				GeneInfoTmp.GeneEnd=End;
				GeneInfoTmp.GeneLength=(End-Start+1);	
				GeneInfoTmp.CDSList.push_back({Start,End});

				if (RefIn)
				{
					for (int ii=Start-1; ii<End ; ii++)
					{
						GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[i])[ii]];
					}
				}

				if (Start==1)
				{
					TmpGene[GeneID]=GeneInfoTmp;
					GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(i,TmpGene));

					map <int,int> AA;
					AA.insert(map <int,int>  :: value_type(Start,End));
					RegionMerger.insert(map <int,map <int,int> >  :: value_type(i,AA));

				}
				else
				{
					GeneDataIT=GeneData.find(i);
					(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;

					MergerIt=RegionMerger.find(i);
					(MergerIt->second).insert(map <int,int>  :: value_type(Start,End)) ;
				}
				End+=2;
			}
		}

	}

	(paraFA04->InStr3)=(paraFA04->InStr3).substr(0,(paraFA04->InStr3).length()-3);
	string path=(paraFA04->InStr3);	
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	string PrefixO=path;

	if (ext == "stat" || ext == "bed")
	{
		PrefixO=path.substr(0,path.rfind('.'));
	}

	string  OutStatFile=PrefixO+".gene.stat.gz";
	string  OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";

	if  ((paraFA04->InInt2)==3)
	{
		OutStatFile=PrefixO+".bed.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tRegionID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==4)
	{
		OutStatFile=PrefixO+".bed.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==5 ||  (paraFA04->InInt2)==6 )
	{
		OutStatFile=PrefixO+".win.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==0)
	{
		OutStatFile=PrefixO+".chr.stat.gz";
		OutHeader="#Chr\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}

	if (RefIn)
	{
		RefBase.clear();
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		if  ((paraFA04->InInt2)==3)
		{
			OutHeader="#Chr\tStart\tEnd\tRegionID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==4)
		{
			OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==5|| (paraFA04->InInt2)==6)
		{
			OutHeader="#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==0)
		{
			OutHeader="#Chr\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
	}
	ogzstream  OUT (OutStatFile.c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<OutStatFile<<endl;
		delete  paraFA04 ; return  0;
	}

	SiteInfo **depth = new SiteInfo *[n_targets];

	for(int i = 0; i < (n_targets); i++)  
	{
		int CC=target_len[i]+500;

		MergerIt=RegionMerger.find(i);
		if (MergerIt==RegionMerger.end())
		{
			CC=500;
		}

		depth[i] = new SiteInfo  [CC];

		for (int32_t j =0 ; j< CC ; j++)
		{
			depth[i][j].Depth=0;
		}
	}



	vector<string> inf; 
	int tid=0;
	int qual=0;
	string line;
	std::string substring = "tp:A:S";

	for (auto it = (paraFA04->InStr1List).begin(); it != (paraFA04->InStr1List).end(); ++it) 
	{
		BamPath=*it;
		igzstream INB (BamPath.c_str(),ifstream::in);
		if(!INB.good())
		{
			cerr << "open PAF file error: "<<BamPath<<endl;
			return 1;
		}

		while(!INB.eof())
		{

			getline(INB,line);
			if (line.length()<=0)  { continue  ; }

			if (  0x100 &  (paraFA04->flags) )
			{
				if (line.find(substring) != std::string::npos)
				{
					continue ;
				}
			}

			inf.clear();
			split(line,inf," \t");
			tid=Chr2IntMap[inf[5]];

			qual=atoi(inf[11].c_str());
			if ( qual < (paraFA04->InInt) )
			{
				continue ;
			}

			int32_t StartRead=atoi(inf[7].c_str()) ;
			int32_t EndRead=atoi(inf[8].c_str()) ;
			if  (StartRead> EndRead)
			{
				int32_t tmp=StartRead;
				StartRead=EndRead;
				EndRead=tmp;
			}

			int indxcg=findFirstStringIndex(inf);
			if (indxcg>1)
			{
				std::vector<char> charVec;
				std::vector<int> intVec;
				parcigar( inf[indxcg] ,  charVec,  intVec);
				int n_cigar=intVec.size();
				int  endTmp ;
				for(int i=0; i <n_cigar;++i)
				{
					int cig=charVec[i];
					int ncig = intVec[i];

					switch (cig)
					{
						case 'M':
						case '=':
						case 'X' :
							endTmp=StartRead+ncig;
							for (  ; StartRead<endTmp;StartRead++)
							{
								(depth[tid][StartRead]).Depth++;
							}
							break;
						case 'D':
						case 'N':
							StartRead=StartRead+ncig;
							break;
					}
				}
			}
			else
			{
				for(int i=StartRead-1; i <EndRead;++i)
				{
					(depth[tid][i]).Depth++;
				}
			}
		}

		INB.close();
	}



	for(int i = 0; i < (n_targets); i++)  
	{
		MergerIt=RegionMerger.find(i);
		if (MergerIt==RegionMerger.end())
		{

		}
		else
		{
			GeneDataIT=GeneData.find(i);
			StatChrDepthLowMEM (depth[i] , GeneDataIT->second, paraFA04->minDep);
		}
	}



	if(paraFA04->SiteOutPut)
	{
		string  OutSSiteFile=PrefixO+".SiteDepth.gz";
		ogzstream  OUTFA (OutSSiteFile.c_str());
		for(int i = 0; i < (n_targets); i++)  
		{
			int CC=(target_len[i]);
			MergerIt=RegionMerger.find(i);
			if (MergerIt==RegionMerger.end())
			{
				CC=100;
				continue ;
			}
			string ChrName=target_name[i];
			for (int32_t j =0 ; j< CC ; j++)
			{
				OUTFA<<ChrName<<"\t"<<j<<"\t"<<depth[i][j].Depth<<"\n";
			}
		}
		OUTFA.close() ;
	}

	if (   (paraFA04->InInt2)==6 )
	{
		ubit64_t  SS_Len =0;
		ubit64_t  SS_Cov=0;
		ubit64_t  SS_TotalD=0;
		if (RefIn)
		{
			OUT<<OutHeader;
			ubit64_t SS_GCGC =0;
			for(int i = 0; i < (n_targets); i++)  
			{
				int CC=(target_len[i]);
				MergerIt=RegionMerger.find(i);
				if (MergerIt==RegionMerger.end())
				{
					CC=100;
					continue ;
				}
				string ChrName=target_name[i];
				int Start ; int End;
				int GeneLength ; int GeneCover ; int GeneDepth=0;
				int GeneGCGC=0;

				for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
				{

					Start=j-1; End=Start+(paraFA04->WinSize);
					if (End>CC ) {End=CC;}
					GeneCover=0;GeneDepth=0;GeneGCGC=0;

					for ( ;Start<End ; Start++)
					{
						if (depth[i][Start].Depth >= paraFA04->minDep)
						{
							GeneCover++;
							GeneDepth+=(depth[i][Start].Depth);
						}
						GeneGCGC+=GCGCArry[(RefBase[i])[Start]];
					}
					GeneLength=End-j+1;
					double Coverage=(GeneCover)*100.0/(GeneLength);
					double MeanDepth=(GeneDepth)*1.0/(GeneLength);
					double GeneGC=(GeneGCGC)*100.0/(GeneLength);
					OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\n";

					SS_Cov+=(GeneCover);
					SS_Len+=(GeneLength);
					SS_GCGC+=(GeneGCGC);
					SS_TotalD+=(GeneDepth);
				}
			}
			double Coverage=SS_Cov*100.0/SS_Len;
			double MeanDepth=SS_TotalD*1.0/SS_Len;
			double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

			OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tGC(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ALLGeneGCRation<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;


		}
		else
		{

			OUT<<OutHeader;
			for(int i = 0; i < (n_targets); i++)  
			{
				int CC=(target_len[i]);
				MergerIt=RegionMerger.find(i);
				if (MergerIt==RegionMerger.end())
				{
					CC=100;
					continue ;
				}
				string ChrName=target_name[i];
				int Start ; int End;
				int GeneLength ; int GeneCover ; int GeneDepth=0;

				for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
				{

					Start=j-1; End=Start+(paraFA04->WinSize);
					if (End>CC ) {End=CC;}
					GeneCover=0;GeneDepth=0;

					for ( ;Start<End ; Start++)
					{
						if (depth[i][Start].Depth >= paraFA04->minDep )
						{
							GeneCover++;
							GeneDepth+=(depth[i][Start].Depth);
						}
					}
					GeneLength=End-j+1;
					double Coverage=(GeneCover)*100.0/(GeneLength);
					double MeanDepth=(GeneDepth)*1.0/(GeneLength);
					OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<"\n";

					SS_Cov+=(GeneCover);
					SS_Len+=(GeneLength);
					SS_TotalD+=(GeneDepth);
				}
			}
			double Coverage=SS_Cov*100.0/SS_Len;
			double MeanDepth=SS_TotalD*1.0/SS_Len;

			OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;
		}
	}

	for(int i = 0; i <(n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;

	cout <<"INFO: Input data read done"<<endl;

	ubit64_t  SS_Len =0;
	ubit64_t  SS_Cov=0;
	ubit64_t  SS_TotalD=0;

	if ((RefIn)  &&  ((paraFA04->InInt2)!=6))
	{
		OUT<<OutHeader;
		ubit64_t  SS_GCGC=0;
		if  ((paraFA04->InInt2)==1 ||  (paraFA04->InInt2)==2   ||  (paraFA04->InInt2)==3   ||  (paraFA04->InInt2)==4 )
		{
			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 
				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);
					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<GeneID<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}
				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		else if ((paraFA04->InInt2)==0)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=target_name[ChrNum];

				ubit64_t  ChrSS_Len =0;
				ubit64_t  ChrSS_Cov=0;
				ubit64_t  ChrSS_TotalD=0;
				ubit64_t  ChrSS_GCGC=0;

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					ChrSS_Cov+=((GeneInfoIT->second).GeneCover);
					ChrSS_Len+=((GeneInfoIT->second).GeneLength);
					ChrSS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					ChrSS_TotalD+=((GeneInfoIT->second).GeneDepth);
				}

				double GeneGC=ChrSS_GCGC*100.0/ChrSS_Len;
				double MeanDepth=ChrSS_TotalD*1.0/ChrSS_Len;
				double Coverage=ChrSS_Cov*100.0/ChrSS_Len;
				OUT<<ChrName<<"\t"<<ChrSS_Len<<"\t"<<ChrSS_Cov<<"\t"<<ChrSS_TotalD<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\n";

			}
		}
		else  if ((paraFA04->InInt2)==5)
		{


			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		double Coverage=SS_Cov*100.0/SS_Len;
		double MeanDepth=SS_TotalD*1.0/SS_Len;
		double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

		OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tGC(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ALLGeneGCRation<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;

	}
	else if (((paraFA04->InInt2)!=6))
	{
		OUT<<OutHeader;
		if  ((paraFA04->InInt2)==1 ||  (paraFA04->InInt2)==2   ||  (paraFA04->InInt2)==3   ||  (paraFA04->InInt2)==4 )
		{
			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);
					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<GeneID<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth;

					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}

				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		else if ((paraFA04->InInt2)==0)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=target_name[ChrNum];

				ubit64_t  ChrSS_Len =0;
				ubit64_t  ChrSS_Cov=0;
				ubit64_t  ChrSS_TotalD=0;

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					ChrSS_Cov+=((GeneInfoIT->second).GeneCover);
					ChrSS_Len+=((GeneInfoIT->second).GeneLength);
					ChrSS_TotalD+=((GeneInfoIT->second).GeneDepth);
				}

				double MeanDepth=ChrSS_TotalD*1.0/ChrSS_Len;
				double Coverage=ChrSS_Cov*100.0/ChrSS_Len;
				OUT<<ChrName<<"\t"<<ChrSS_Len<<"\t"<<ChrSS_Cov<<"\t"<<ChrSS_TotalD<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<"\n";

			}
		}
		else if ((paraFA04->InInt2)==5)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		double Coverage=SS_Cov*100.0/SS_Len;
		double MeanDepth=SS_TotalD*1.0/SS_Len;

		OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tCoverage(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;

	}

	OUT.close();

	return 0;
}




int BamList_main(In3str1v *paraFA04 )
{


	string  BamPath=(paraFA04->InStr1);
	if (BamPath.length()<=3)
	{
		cerr<<"Error: Failed to open the bam file: "<<BamPath<<endl;
		return 1;
	}

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	hts_set_log_level(HTS_LOG_OFF);
	//
	if(BamIn->format.format == htsExactFormat::cram)
	{
		const char* ref_file = (paraFA04->reference).c_str();
		hts_set_fai_filename(BamIn, ref_file);
		hts_set_opt(BamIn, CRAM_OPT_DECODE_MD, 0);
		hts_set_opt(BamIn, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	}
	//
	int numThreads = (paraFA04->CPU);
	hts_set_opt(BamIn, HTS_OPT_NTHREADS, numThreads);
	header = sam_hdr_read(BamIn);

	map <string,int> Chr2IntMap; 
	for(int i = 0; i < (header->n_targets); i++) 
	{
		string ChrName=header->target_name[i];
		Chr2IntMap.insert( map <string,int>  :: value_type (ChrName,i));
	}

	bam1_t *alnTA = bam_init1();

	sam_close(BamIn);

	map <int,string>  RefBase ;
	bool  RefIn=false;
	if ((paraFA04->gc)==true)
	{
		if (!(paraFA04->reference).empty()){
			gzFile fpRef;
			kseq_t *seq;
			int l;
			fpRef = gzopen((paraFA04->reference).c_str(), "r");
			seq = kseq_init(fpRef);
			RefIn=true;
			while ((l = kseq_read(seq)) >= 0 )
			{
				string chr=(seq->name.s);
				int ID=Chr2IntMap[chr];
				string seqBB=seq->seq.s;
				RefBase.insert( map <int,string>  :: value_type (ID,seqBB));
			}
			kseq_destroy(seq); gzclose(fpRef);
		}
		else{
			cerr<< "Error: lack reference sequence (-r) for GC parse"<<endl ;
			return 0;
		}
	}
	//
	int GCGCArry[256] = {0};
	GCGCArry['C']=1; GCGCArry['c']=1; GCGCArry['G']=1; GCGCArry['g']=1; 
	GCGCArry['a']=0; GCGCArry['A']=0; GCGCArry['t']=0; GCGCArry['T']=0;

	map <string,int> :: iterator  MapItChr2Int ;
	unordered_map <string,int> :: iterator  UnMapItChr2Int ;
	int Start ; int End ;  string ChrName ;
	string GeneID ;
	map <int,map <string,GeneInfo> > GeneData;
	map <int,map <string,GeneInfo> >  :: iterator  GeneDataIT;
	map <string,GeneInfo>   :: iterator  GeneInfoIT  ;

	if  ((paraFA04->InInt2)!=0)
	{
		igzstream LIST ((paraFA04->InStr2).c_str(),ifstream::in); // igzstream
		if (!LIST.good())
		{
			cerr << "Error: Cannot open the GFF/GTF File: "<<(paraFA04->InStr1)<<endl;
			return  1;
		}

		if ((paraFA04->InInt2)==1)
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				istringstream isone (line,istringstream::in);
				string flag , CDS , ZhengFu ,geneID ;;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End >>flag>>ZhengFu>>flag>>geneID;;

				vector<string> inf;
				vector<string> Temp;
				split(geneID,inf,",;");
				split(inf[0],Temp,"=");
				string GeneID= Temp[Temp.size()-1];
				for ( int jk=1;  jk<inf.size() ; jk++)
				{
					vector<string> Temptmp2;
					split(inf[jk],Temptmp2,"=");
					if (Temptmp2[0] == "Parent")
					{
						GeneID= Temptmp2[Temptmp2.size()-1];
					}
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{
					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];
							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}

				}
			}

		}
		else if ((paraFA04->InInt2)==2)
		{
			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				line=replace_all(line,"\"","");
				line=replace_all(line,";","");

				istringstream isone (line,istringstream::in);
				string  flag , CDS , ZhengFu ,geneID;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End >>flag>>ZhengFu>>flag  ;

				vector<string> inf;
				split(line,inf,"\t ");
				string GeneID=inf[9];
				for ( int jk=8;  jk<inf.size() ; jk+=2)
				{

					if (inf[0] == "transcript_id")
					{
						GeneID= inf[jk+1];
					}
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{
					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}


				}
			}
		}
		else if ((paraFA04->InInt2)==3)
		{
			string StartStr , EndStr ;
			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);
				if (line.length()<=0)  {continue;}
				if (line[0] == '#')  { continue;}
				istringstream isone (line,istringstream::in);
				isone>>ChrName>>StartStr>>EndStr;
				GeneID=ChrName+"_"+StartStr+"_"+EndStr ;
				Start=atoi(StartStr.c_str()); End=atoi(EndStr.c_str());
				if  (Start > End )
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
					continue;
				}
				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{

					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}
				}
			} // LIST.eof
		}  // end (paraFA04->InInt2)==3

		else if ((paraFA04->InInt2)==4)
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);
				if (line.length()<=0)  {continue;}
				if (line[0] == '#')  { continue;}
				istringstream isone (line,istringstream::in);
				isone>>ChrName>>Start>>End>>GeneID;
				if  (Start > End )
				{
					cerr<<line<<"Warning: This region may be incorrect. \n"<<endl;
					continue;
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{

					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn)
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}
				}
			}
		}
		else
		{

		}

		LIST.close();
	}



	map <int,map <int,int> > RegionMerger;

	map <int,map <int,int> > :: iterator MergerIt ;
	map <int,int> :: iterator  MergerMapSSEE  ;



	for ( GeneDataIT = GeneData.begin(); GeneDataIT != GeneData.end(); ++GeneDataIT)
	{

		int ChrInt = GeneDataIT ->first ;
		map <int,int>  RegionAAA;
		map <int,int> :: iterator MapSSEE ;
		for (GeneInfoIT=(GeneDataIT -> second).begin() ; GeneInfoIT!=(GeneDataIT -> second).end()  ;  GeneInfoIT ++)
		{
			Start=(GeneInfoIT->second).GeneStart;
			End=(GeneInfoIT->second).GeneEnd;
			MapSSEE=RegionAAA.find(Start);
			if (MapSSEE==RegionAAA.end())
			{
				RegionAAA.insert(map <int,int>  :: value_type(Start,End)) ;
			}
			else
			{
				if (End >(MapSSEE->second))
				{
					MapSSEE->second=End;
				}
			}
		}




		MapSSEE=RegionAAA.begin() ;
		Start=MapSSEE->first;
		End=MapSSEE->second;

		map <int,int> Start2End ;
		Start2End[Start]=End;
		RegionMerger.insert(map <int,map <int,int> > ::value_type(ChrInt,Start2End));
		MergerIt=RegionMerger.find(ChrInt);
		MapSSEE++;


		for(  ; MapSSEE!=RegionAAA.end() ; MapSSEE++ )
		{
			if (  (MapSSEE->first) >  End    )
			{
				Start=MapSSEE->first;
				End=MapSSEE->second;
				(MergerIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
			}
			else if  (  (MapSSEE->second) >  End    )
			{
				MergerMapSSEE=(MergerIt->second).end(); MergerMapSSEE--;
				MergerMapSSEE->second= (MapSSEE->second) ;
				End=MapSSEE->second;
			}
		}
	}







	// check region
	if (RegionMerger.empty())
	{
		//MeMBinWindows
		int  MeMBinWindows=10000000;
		if ((paraFA04->WinSize) ==0 )
		{
			//cout<<"Warning: GFF/GTF or BED was not provided, the total chromosome length will be used as the parsing region.\n";
			(paraFA04->InInt2)=0;
		}
		else if ((paraFA04->WinSize) <150 )
		{
			(paraFA04->InInt2)=6;
		}
		else
		{
			(paraFA04->InInt2)=5;
			MeMBinWindows=(paraFA04->WinSize);
		}


		for(int i = 0; i < (header->n_targets); i++)
		{
			Start=1;
			End=2;

			for (Start=1;  End <= (header->target_len[i]) ; Start+=MeMBinWindows )
			{
				End=Start+MeMBinWindows-1;
				if (End>(header->target_len[i]))
				{
					End=(header->target_len[i]);
				}


				string GeneID=header->target_name[i]+Int2Str(Start);
				map <string,GeneInfo> TmpGene;
				GeneInfo  GeneInfoTmp ;
				GeneInfoTmp.GeneStart=Start;
				GeneInfoTmp.GeneEnd=End;
				GeneInfoTmp.GeneLength=(End-Start+1);	
				GeneInfoTmp.CDSList.push_back({Start,End});

				if (RefIn)
				{
					for (int ii=Start-1; ii<End ; ii++)
					{
						GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[i])[ii]];
					}
				}

				if (Start==1)
				{
					TmpGene[GeneID]=GeneInfoTmp;
					GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(i,TmpGene));

					map <int,int> AA;
					AA.insert(map <int,int>  :: value_type(Start,End));
					RegionMerger.insert(map <int,map <int,int> >  :: value_type(i,AA));

				}
				else
				{
					GeneDataIT=GeneData.find(i);
					(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;

					MergerIt=RegionMerger.find(i);
					(MergerIt->second).insert(map <int,int>  :: value_type(Start,End)) ;
				}

				End+=2;

			}


		}

	}





	(paraFA04->InStr3)=(paraFA04->InStr3).substr(0,(paraFA04->InStr3).length()-3);
	string path=(paraFA04->InStr3);	
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	string PrefixO=path;

	if (ext == "stat" || ext == "bed")
	{
		PrefixO=path.substr(0,path.rfind('.'));
	}

	string  OutStatFile=PrefixO+".gene.stat.gz";
	string  OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";

	if  ((paraFA04->InInt2)==3)
	{
		OutStatFile=PrefixO+".bed.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tRegionID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==4)
	{
		OutStatFile=PrefixO+".bed.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==5 ||  (paraFA04->InInt2)==6 )
	{
		OutStatFile=PrefixO+".win.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==0)
	{
		OutStatFile=PrefixO+".chr.stat.gz";
		OutHeader="#Chr\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}

	if (RefIn)
	{
		RefBase.clear();
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		if  ((paraFA04->InInt2)==3)
		{
			OutHeader="#Chr\tStart\tEnd\tRegionID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==4)
		{
			OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==5|| (paraFA04->InInt2)==6)
		{
			OutHeader="#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==0)
		{
			OutHeader="#Chr\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
	}
	ogzstream  OUT (OutStatFile.c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<OutStatFile<<endl;
		return  0;
	}







	SiteInfo **depth = new SiteInfo *[(header->n_targets)];
	uint32_t flags = (paraFA04->flags);

	for(int i = 0; i < (header->n_targets); i++)  
	{
		int CC=(header->target_len[i])+100;

		depth[i] = new SiteInfo [CC];
		for (int32_t j =0 ; j< CC ; j++)
		{
			depth[i][j].Depth=0;
		}
	}




	for (auto it = (paraFA04->InStr1List).begin(); it != (paraFA04->InStr1List).end(); ++it) 
	{
		BamPath=*it;
		// read index file
		string bambai=BamPath+".bai";
		string bamcsi=BamPath+".csi";
		string crambai=BamPath+".crai";


		if ((( access(bambai.c_str(), 0) == 0 )  ||  (access(crambai.c_str(), 0) == 0 )   ||   (access(bamcsi.c_str(), 0) == 0 )   )  && (paraFA04->TF ) )	
		{

			ubit64_t GenomeLen=0;
			vector< pair<int, int> > Int2Len;
			auto RegionIt=RegionMerger.begin();
			int ChrNumRun=0;
			int MaxChrLenTen=0;
			while (RegionIt!= RegionMerger.end())
			{
				int  ChrNum=RegionIt->first;
				GenomeLen+=(header->target_len[ChrNum]);
				Int2Len.push_back({ChrNum,header->target_len[ChrNum]});
				if (header->target_len[ChrNum]>MaxChrLenTen)
				{
					MaxChrLenTen=header->target_len[ChrNum];
				}
				RegionIt++;
				ChrNumRun++;
			}


			MaxChrLenTen=int(MaxChrLenTen/10);
			RegionIt=RegionMerger.begin();
			int BigChrNum=0;
			while (RegionIt!= RegionMerger.end())
			{
				int  ChrNum=RegionIt->first;
				if (header->target_len[ChrNum]>MaxChrLenTen)
				{
					BigChrNum++;
				}
				RegionIt++;
			}

			int BamThreadSum=0;
			if ( (paraFA04->CPU) >  BigChrNum )
			{
				BamThreadSum=(paraFA04->CPU)-BigChrNum;
				(paraFA04->CPU)=BigChrNum;
			}
			ubit64_t meanLenDea=GenomeLen/(paraFA04->CPU)+2;
			map <int,vector <int> >  VecChr;
			map <int,vector <int> > :: iterator  VecChrITT;

			sort(Int2Len.begin(),Int2Len.end(),[](const pair<int, int>& a, const pair<int, int>& b){ return a.second > b.second; }); 

			int StBase=0;
			int ShiftQ=0;
			int RunCPU=(paraFA04->CPU);
			ubit64_t *CutThread = new ubit64_t [RunCPU];
			int  *BamThread =new int[RunCPU];
			for (int i = 0; i < RunCPU ;++i)
			{
				CutThread[i]=0;
				BamThread[i]=0;
			}

			while(BamThreadSum>0)
			{
				ShiftQ=StBase%RunCPU;
				StBase++;
				BamThread[ShiftQ]++;
				BamThreadSum--;
			}
			ShiftQ=0;
			StBase=0;
			for (int i = 0; i != Int2Len.size(); i++)
			{
				int  ChrNum=Int2Len[i].first;
				int  Num=int((i-ShiftQ)%RunCPU)+StBase;

				VecChrITT=VecChr.find(Num);
				if  (VecChrITT==VecChr.end())
				{
					vector <int> tmp;
					tmp.emplace_back(ChrNum);
					VecChr.insert( map <int,vector <int> >  :: value_type (Num,tmp));
				}
				else
				{
					(VecChrITT->second).emplace_back(ChrNum);
				}

				CutThread[Num]+=Int2Len[i].second;
				if (  CutThread[Num] >meanLenDea    )
				{
					StBase++;
					RunCPU--;
					ShiftQ=i+1;
				}
			}
			delete [] CutThread ;

			vector<thread> ThreadsVector ;
			VecChrITT=VecChr.begin();
			ShiftQ=0;
			while(VecChrITT!=VecChr.end())
			{
				if  (BamThread[ShiftQ]> 0 ) {BamThread[ShiftQ]=2;}
				ThreadsVector.emplace_back(ProDealChrBambaiALLSite , std::ref(BamPath), paraFA04,std::ref(depth) ,std::ref(RegionMerger) , std::ref(VecChrITT->second),std::ref(GeneData) ,std::ref(BamThread[ShiftQ]));
				VecChrITT++;
				ShiftQ++;
			}
			int AA=ThreadsVector.size();
			for (int i=0 ; i<AA ; i++)
			{
				ThreadsVector[i].join();
			}
			ThreadsVector.clear();
			delete [] BamThread ;
		}
		else
		{
			bam1_t *aln = bam_init1();
			uint32_t *cigar;

			bam_hdr_t *headerRR;
			samFile *BamInRR = hts_open(BamPath.c_str(), "r");
			int numThreads = (paraFA04->CPU);
			hts_set_opt(BamInRR, HTS_OPT_NTHREADS, numThreads);
			if(BamInRR->format.format == htsExactFormat::cram)
			{
				const char* ref_file = (paraFA04->reference).c_str();
				hts_set_fai_filename(BamInRR, ref_file);
				hts_set_opt(BamInRR, CRAM_OPT_DECODE_MD, 0);
				hts_set_opt(BamInRR, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ| SAM_CIGAR);
			}
			headerRR = sam_hdr_read(BamInRR);


			string headertext=headerRR->text;
			string::size_type  posA = (headertext).find("\tSO:");
			bool sortBam=false;

			if (posA!=string::npos)
			{
				posA+=4;
				string::size_type  posB=(headertext).find_first_of("\n\t", posA);
				string  sortinfo=(headertext).substr(posA, posB-posA);
				if (sortinfo=="coordinate")
				{
					sortBam=true;
				}
			}

			if (sortBam)
			{
				cout <<"Warning: PanDepth will run in No Index mode: "<<BamPath<<endl;
				bool ALLEnd=false;
				bool  *EndChr = new bool [(header->n_targets)];
				map <int,int> :: iterator  *ArryIt = new map <int,int> :: iterator [(header->n_targets)];
				map <int,int> :: iterator  *ArryItEnd = new map <int,int> :: iterator [(header->n_targets)];

				for(int i = 0; i < (header->n_targets); i++)
				{
					EndChr[i]=false ;
					MergerIt=RegionMerger.find(i);
					if (MergerIt==RegionMerger.end())
					{
						EndChr[i]=true ;
					}
					else
					{
						ArryIt[i]=(MergerIt->second).begin();
						ArryItEnd[i]=(MergerIt->second).end();
					}
				}

				while (sam_read1(BamInRR, headerRR, aln) >= 0)
				{
					if ( EndChr[(aln->core).tid])  {continue ;}	
					if ( (aln->core).qual < (paraFA04->InInt) )
					{
						continue ;
					}
					if ( aln->core.flag & flags ) continue;
					int32_t EndRead=bam_endpos(aln);
					if (EndRead <(ArryIt[(aln->core).tid]->first)) {continue;}

					int32_t StartRead=((aln->core).pos);
					if ( StartRead  > (ArryIt[(aln->core).tid]->second) )
					{

						for (ArryIt[(aln->core).tid]++  ; ArryIt[(aln->core).tid]!=ArryItEnd[(aln->core).tid] ;ArryIt[(aln->core).tid]++)
						{
							if ( StartRead <=  (ArryIt[(aln->core).tid]->second)  )
							{
								break;
							}
						}
						if  (ArryIt[(aln->core).tid]==ArryItEnd[(aln->core).tid])
						{
							EndChr[(aln->core).tid]=true;
							ALLEnd=true;
							for(int i = 0; i < (header->n_targets); i++)
							{
								if (!EndChr[i])
								{
									ALLEnd=false;
								}
							}
							if (ALLEnd)
							{
								break;
							}
						}
					}


					cigar = bam_get_cigar(aln);
					int  endTmp ;
					for(int i=0; i < aln->core.n_cigar;++i)
					{
						int cig=bam_cigar_op(cigar[i]);
						int ncig = bam_cigar_oplen(cigar[i]);
						switch (cig)
						{
							case 0:
							case 7:
							case 8 :
								endTmp=StartRead+ncig;
								for (  ; StartRead<endTmp;StartRead++)
								{
									(depth[(aln->core).tid][StartRead]).Depth++;
								}
								break;
							case 2:
							case 3:
								StartRead=StartRead+ncig;
								break;
						}
					}
				}
				delete[] EndChr;
				delete[] ArryItEnd;
			}
			else
			{
				cout <<"Warning: Can't find index file of input BAM/CRAM. PanDepth will run in No Index mode: "<<BamPath<<endl;
				while (sam_read1(BamInRR, headerRR, aln) >= 0)
				{

					if ( (aln->core).qual < (paraFA04->InInt))
					{
						continue ;
					}
					if ( aln->core.flag & flags ) continue;
					int32_t EndRead=bam_endpos(aln);
					int32_t StartRead=((aln->core).pos);
					cigar = bam_get_cigar(aln);
					int  endTmp ;
					for(int i=0; i < aln->core.n_cigar;++i)
					{
						int cig=bam_cigar_op(cigar[i]);
						int ncig = bam_cigar_oplen(cigar[i]);
						switch (cig)
						{
							case 0:
							case 7:
							case 8 :
								endTmp=StartRead+ncig;
								for (  ; StartRead<endTmp;StartRead++)
								{
									(depth[(aln->core).tid][StartRead]).Depth++;
								}
								break;
							case 2:
							case 3:
								StartRead=StartRead+ncig;
								break;
						}
					}
				}
			}

			sam_close(BamInRR);
			bam_hdr_destroy(headerRR);
			bam_destroy1(aln);


		}



	}



	for(int i = 0; i < (header->n_targets); i++)  
	{
		MergerIt=RegionMerger.find(i);
		if (MergerIt==RegionMerger.end())
		{

		}
		else
		{
			GeneDataIT=GeneData.find(i);
			StatChrDepthLowMEM ( depth[i] , GeneDataIT->second, paraFA04->minDep);
		}
	}


	if(paraFA04->SiteOutPut)
	{
		string  OutSSiteFile=PrefixO+".SiteDepth.gz";
		ogzstream  OUTFA (OutSSiteFile.c_str());
		for(int i = 0; i < (header->n_targets); i++)  
		{
			int CC=(header->target_len[i]);
			MergerIt=RegionMerger.find(i);
			if (MergerIt==RegionMerger.end())
			{
				CC=100;
				continue ;
			}
			string ChrName=header->target_name[i];
			for (int32_t j =0 ; j< CC ; j++)
			{
				OUTFA<<ChrName<<"\t"<<j<<"\t"<<depth[i][j].Depth<<"\n";
			}
		}
		OUTFA.close() ;
	}


	if (   (paraFA04->InInt2)==6 )
	{
		ubit64_t  SS_Len =0;
		ubit64_t  SS_Cov=0;
		ubit64_t  SS_TotalD=0;
		if (RefIn)
		{
			OUT<<OutHeader;
			ubit64_t SS_GCGC =0;
			for(int i = 0; i < (header->n_targets); i++)  
			{
				int CC=(header->target_len[i]);
				MergerIt=RegionMerger.find(i);
				if (MergerIt==RegionMerger.end())
				{
					CC=100;
					continue ;
				}
				string ChrName=header->target_name[i];
				int Start ; int End;
				int GeneLength ; int GeneCover ; int GeneDepth=0;
				int GeneGCGC=0;

				for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
				{

					Start=j-1; End=Start+(paraFA04->WinSize);
					if (End>CC ) {End=CC;}
					GeneCover=0;GeneDepth=0;GeneGCGC=0;

					for ( ;Start<End ; Start++)
					{
						if (depth[i][Start].Depth >= paraFA04->minDep )
						{
							GeneCover++;
							GeneDepth+=(depth[i][Start].Depth);
						}
						GeneGCGC+=GCGCArry[(RefBase[i])[Start]];
					}
					GeneLength=End-j+1;
					double Coverage=(GeneCover)*100.0/(GeneLength);
					double MeanDepth=(GeneDepth)*1.0/(GeneLength);
					double GeneGC=(GeneGCGC)*100.0/(GeneLength);
					OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\n";

					SS_Cov+=(GeneCover);
					SS_Len+=(GeneLength);
					SS_GCGC+=(GeneGCGC);
					SS_TotalD+=(GeneDepth);
				}
			}
			double Coverage=SS_Cov*100.0/SS_Len;
			double MeanDepth=SS_TotalD*1.0/SS_Len;
			double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

			OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tGC(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ALLGeneGCRation<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;


		}
		else
		{

			OUT<<OutHeader;
			for(int i = 0; i < (header->n_targets); i++)  
			{
				int CC=(header->target_len[i]);
				MergerIt=RegionMerger.find(i);
				if (MergerIt==RegionMerger.end())
				{
					CC=100;
					continue ;
				}
				string ChrName=header->target_name[i];
				int Start ; int End;
				int GeneLength ; int GeneCover ; int GeneDepth=0;

				for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
				{

					Start=j-1; End=Start+(paraFA04->WinSize);
					if (End>CC ) {End=CC;}
					GeneCover=0;GeneDepth=0;

					for ( ;Start<End ; Start++)
					{
						if (depth[i][Start].Depth >= paraFA04->minDep )
						{
							GeneCover++;
							GeneDepth+=(depth[i][Start].Depth);
						}
					}
					GeneLength=End-j+1;
					double Coverage=(GeneCover)*100.0/(GeneLength);
					double MeanDepth=(GeneDepth)*1.0/(GeneLength);
					OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<"\n";

					SS_Cov+=(GeneCover);
					SS_Len+=(GeneLength);
					SS_TotalD+=(GeneDepth);
				}
			}
			double Coverage=SS_Cov*100.0/SS_Len;
			double MeanDepth=SS_TotalD*1.0/SS_Len;

			OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;
		}
	}


	for(int i = 0; i <(header->n_targets); i++)
	{
		delete[] depth[i];  
	}
	delete[] depth;


	cout <<"INFO: Input data read done"<<endl;

	ubit64_t  SS_Len =0;
	ubit64_t  SS_Cov=0;
	ubit64_t  SS_TotalD=0;

	if ((RefIn)  &&  ((paraFA04->InInt2)!=6))
	{
		OUT<<OutHeader;
		ubit64_t  SS_GCGC=0;
		if  ((paraFA04->InInt2)==1 ||  (paraFA04->InInt2)==2   ||  (paraFA04->InInt2)==3   ||  (paraFA04->InInt2)==4 )
		{
			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 
				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);
					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<GeneID<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}
				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		else if ((paraFA04->InInt2)==0)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];

				ubit64_t  ChrSS_Len =0;
				ubit64_t  ChrSS_Cov=0;
				ubit64_t  ChrSS_TotalD=0;
				ubit64_t  ChrSS_GCGC=0;

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					ChrSS_Cov+=((GeneInfoIT->second).GeneCover);
					ChrSS_Len+=((GeneInfoIT->second).GeneLength);
					ChrSS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					ChrSS_TotalD+=((GeneInfoIT->second).GeneDepth);
				}

				double GeneGC=ChrSS_GCGC*100.0/ChrSS_Len;
				double MeanDepth=ChrSS_TotalD*1.0/ChrSS_Len;
				double Coverage=ChrSS_Cov*100.0/ChrSS_Len;
				OUT<<ChrName<<"\t"<<ChrSS_Len<<"\t"<<ChrSS_Cov<<"\t"<<ChrSS_TotalD<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\n";

			}
		}
		else  if ((paraFA04->InInt2)==5)
		{


			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		double Coverage=SS_Cov*100.0/SS_Len;
		double MeanDepth=SS_TotalD*1.0/SS_Len;
		double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

		OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tGC(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ALLGeneGCRation<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;

	}
	else if (((paraFA04->InInt2)!=6))
	{
		OUT<<OutHeader;
		if  ((paraFA04->InInt2)==1 ||  (paraFA04->InInt2)==2   ||  (paraFA04->InInt2)==3   ||  (paraFA04->InInt2)==4 )
		{
			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);
					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<GeneID<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth;

					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}

				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		else if ((paraFA04->InInt2)==0)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];

				ubit64_t  ChrSS_Len =0;
				ubit64_t  ChrSS_Cov=0;
				ubit64_t  ChrSS_TotalD=0;

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					ChrSS_Cov+=((GeneInfoIT->second).GeneCover);
					ChrSS_Len+=((GeneInfoIT->second).GeneLength);
					ChrSS_TotalD+=((GeneInfoIT->second).GeneDepth);
				}

				double MeanDepth=ChrSS_TotalD*1.0/ChrSS_Len;
				double Coverage=ChrSS_Cov*100.0/ChrSS_Len;
				OUT<<ChrName<<"\t"<<ChrSS_Len<<"\t"<<ChrSS_Cov<<"\t"<<ChrSS_TotalD<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<"\n";

			}
		}
		else if ((paraFA04->InInt2)==5)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		double Coverage=SS_Cov*100.0/SS_Len;
		double MeanDepth=SS_TotalD*1.0/SS_Len;

		OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tCoverage(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;

	}


	OUT.close();


	bam_hdr_destroy(header);

	return 0;
}


int List_main( In3str1v *paraFA04  ) 
{
	string  BamPath=(paraFA04->InStr1);
	if (BamPath.length()<=0)  
	{
		cerr<<"Error: Failed to open the BAM/CRAM file: "<<BamPath<<endl;
		return 1;
	}

	string ext =BamPath.substr(BamPath.rfind('.') ==string::npos ? BamPath.length() : BamPath.rfind('.') + 1);
	if  (ext == "gz")
	{
		string BamPathA=BamPath.substr(0,BamPath.rfind('.'));
		ext =BamPathA.substr(BamPathA.rfind('.') ==string::npos ? BamPathA.length() : BamPathA.rfind('.') + 1);
	}
	if ( (ext == "paf" ) ||  (ext == "PAF" ) )
	{
		cout <<"INFO: Run PAF format data "<<endl;
		paf_main( paraFA04 ) ;
	}
	else
	{
		BamList_main( paraFA04 ) ;
	}
	return 0;
}


int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=-1;
	int FileNum=(bamCov_help01(argc, argv, paraFA04));

	if ((FileNum==0))
	{
		delete paraFA04 ;
		return 0 ;
	}
	else if (FileNum > 1)
	{
		cout <<"INFO: Run multi-file data "<<endl;
		List_main( paraFA04 ) ;
		delete paraFA04 ;
		return 0 ;
	}

	string  BamPath=(paraFA04->InStr1);
	if (BamPath.length()<=0)  
	{
		cerr<<"Error: Failed to open the BAM/CRAM file: "<<BamPath<<endl;
		return 1;
	}

	string ext =BamPath.substr(BamPath.rfind('.') ==string::npos ? BamPath.length() : BamPath.rfind('.') + 1);
	if  (ext == "gz")
	{
		string BamPathA=BamPath.substr(0,BamPath.rfind('.'));
		ext =BamPathA.substr(BamPathA.rfind('.') ==string::npos ? BamPathA.length() : BamPathA.rfind('.') + 1);
	}

	if ( (ext == "paf" ) ||  (ext == "PAF" ) )
	{
		cout <<"INFO: Run paf Format data "<<endl;
		paf_main( paraFA04 ) ;
		delete paraFA04 ;
		return 0;
	}

	bam_hdr_t *header;
	samFile *BamIn = hts_open(BamPath.c_str(), "r");
	hts_set_log_level(HTS_LOG_OFF);
	//
	if(BamIn->format.format == htsExactFormat::cram)
	{
		const char* ref_file = (paraFA04->reference).c_str();
		hts_set_fai_filename(BamIn, ref_file);
		hts_set_opt(BamIn, CRAM_OPT_DECODE_MD, 0);
		hts_set_opt(BamIn, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);
	}
	//
	int numThreads = (paraFA04->CPU);
	hts_set_opt(BamIn, HTS_OPT_NTHREADS, numThreads);
	header = sam_hdr_read(BamIn);

	map <string,int> Chr2IntMap; 
	for(int i = 0; i < (header->n_targets); i++) 
	{
		string ChrName=header->target_name[i];
		Chr2IntMap.insert( map <string,int>  :: value_type (ChrName,i));
	}

	bam1_t *alnTA = bam_init1();

	sam_close(BamIn);

	map <int,string>  RefBase ;
	bool  RefIn=false;
	if ((paraFA04->gc)==true)
	{
		if (!(paraFA04->reference).empty()){
			gzFile fpRef;
			kseq_t *seq;
			int l;
			fpRef = gzopen((paraFA04->reference).c_str(), "r");
			seq = kseq_init(fpRef);
			RefIn=true;
			while ((l = kseq_read(seq)) >= 0 )
			{
				string chr=(seq->name.s);
				int ID=Chr2IntMap[chr];
				string seqBB=seq->seq.s;
				RefBase.insert( map <int,string>  :: value_type (ID,seqBB));
			}
			kseq_destroy(seq); gzclose(fpRef);
		}
		else{
			cerr<< "Error: lack reference sequence (-r) for GC parse"<<endl ;
			return 0;
		}
	}
	//
	int GCGCArry[256] = {0};
	GCGCArry['C']=1; GCGCArry['c']=1; GCGCArry['G']=1; GCGCArry['g']=1; 
	GCGCArry['a']=0; GCGCArry['A']=0; GCGCArry['t']=0; GCGCArry['T']=0;

	map <string,int> :: iterator  MapItChr2Int ;
	unordered_map <string,int> :: iterator  UnMapItChr2Int ;
	int Start ; int End ;  string ChrName ;
	string GeneID ;
	map <int,map <string,GeneInfo> > GeneData;
	map <int,map <string,GeneInfo> >  :: iterator  GeneDataIT;
	map <string,GeneInfo>   :: iterator  GeneInfoIT  ;

	if  ((paraFA04->InInt2)!=0)
	{

		igzstream LIST ((paraFA04->InStr2).c_str(),ifstream::in); // igzstream
		if (!LIST.good())
		{
			cerr << "Error: Cannot open the GFF/GTF File: "<<(paraFA04->InStr1)<<endl;
			return  1;
		}

		if ((paraFA04->InInt2)==1)
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				istringstream isone (line,istringstream::in);
				string flag , CDS , ZhengFu ,geneID ;;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End >>flag>>ZhengFu>>flag>>geneID;;

				vector<string> inf;
				vector<string> Temp;
				split(geneID,inf,",;");
				split(inf[0],Temp,"=");
				string GeneID= Temp[Temp.size()-1];
				for ( int jk=1;  jk<inf.size() ; jk++)
				{
					vector<string> Temptmp2;
					split(inf[jk],Temptmp2,"=");
					if (Temptmp2[0] == "Parent")
					{
						GeneID= Temptmp2[Temptmp2.size()-1];
					}
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{
					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];
							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}

				}
			}

		}
		else if ((paraFA04->InInt2)==2)
		{
			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);

				if (line.length()<=0 || line[0] == '#' )  { continue  ; }
				line=replace_all(line,"\"","");
				line=replace_all(line,";","");

				istringstream isone (line,istringstream::in);
				string  flag , CDS , ZhengFu ,geneID;
				llong Start,End ;
				isone>>ChrName>>flag>>CDS ;
				if (CDS  != (paraFA04->CDS) )  { continue  ; }
				isone>>Start>>End >>flag>>ZhengFu>>flag  ;

				vector<string> inf;
				split(line,inf,"\t ");
				string GeneID=inf[9];
				for ( int jk=8;  jk<inf.size() ; jk+=2)
				{

					if (inf[0] == "transcript_id")
					{
						GeneID= inf[jk+1];
					}
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{
					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}


				}
			}
		}
		else if ((paraFA04->InInt2)==3)
		{
			string StartStr , EndStr ;
			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);
				if (line.length()<=0)  {continue;}
				if (line[0] == '#')  { continue;}
				istringstream isone (line,istringstream::in);
				isone>>ChrName>>StartStr>>EndStr;
				GeneID=ChrName+"_"+StartStr+"_"+EndStr ;
				Start=atoi(StartStr.c_str()); End=atoi(EndStr.c_str());
				if  (Start > End )
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
					continue;
				}
				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{

					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn )
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}
				}
			} // LIST.eof
		}  // end (paraFA04->InInt2)==3

		else if ((paraFA04->InInt2)==4)
		{

			while(!LIST.eof())
			{
				string  line;
				getline(LIST,line);
				if (line.length()<=0)  {continue;}
				if (line[0] == '#')  { continue;}
				istringstream isone (line,istringstream::in);
				isone>>ChrName>>Start>>End>>GeneID;
				if  (Start > End )
				{
					cerr<<line<<"Warning: This region may be incorrect. \n"<<endl;
					continue;
				}

				MapItChr2Int=Chr2IntMap.find(ChrName);

				if (MapItChr2Int==Chr2IntMap.end())
				{
					cerr<<line<<"Warning: This region may be incorrect.\n"<<endl;
				}
				else
				{

					GeneDataIT=GeneData.find(MapItChr2Int->second);
					if (GeneDataIT==GeneData.end())
					{
						map <string,GeneInfo>   TmpGene ;
						GeneInfo  GeneInfoTmp ;
						GeneInfoTmp.GeneStart=Start;
						GeneInfoTmp.GeneEnd=End;
						GeneInfoTmp.GeneLength=End-Start+1;	
						GeneInfoTmp.CDSList.push_back({Start,End});
						if (RefIn )
						{
							for (int ii=Start-1; ii<End ; ii++)
							{
								GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

							}
						}
						TmpGene[GeneID]=GeneInfoTmp;
						GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(MapItChr2Int->second,TmpGene));
					}
					else
					{
						GeneInfoIT=(GeneDataIT->second).find(GeneID);
						if (GeneInfoIT==(GeneDataIT->second).end())
						{
							GeneInfo  GeneInfoTmp ;
							GeneInfoTmp.GeneStart=Start;
							GeneInfoTmp.GeneEnd=End;
							GeneInfoTmp.GeneLength=End-Start+1;	
							GeneInfoTmp.CDSList.push_back({Start,End});
							if (RefIn)
							{
								for (int ii=Start-1; ii<End ; ii++)
								{
									GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[MapItChr2Int->second])[ii]];

								}
							}

							(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;
						}
						else
						{
							if ((GeneInfoIT->second).GeneStart>Start) {(GeneInfoIT->second).GeneStart=Start;}
							if ((GeneInfoIT->second).GeneEnd<End) {(GeneInfoIT->second).GeneEnd=End;}
							((GeneInfoIT->second).GeneLength)+=(End-Start+1);
							(GeneInfoIT->second).CDSList.push_back({Start,End});
						}
					}
				}
			}
		}
		else
		{

		}

		LIST.close();
	}






	map <int,map <int,int> > RegionMerger;

	map <int,map <int,int> > :: iterator MergerIt ;
	map <int,int> :: iterator  MergerMapSSEE  ;



	for ( GeneDataIT = GeneData.begin(); GeneDataIT != GeneData.end(); ++GeneDataIT)
	{

		int ChrInt = GeneDataIT ->first ;
		map <int,int>  RegionAAA;
		map <int,int> :: iterator MapSSEE ;
		for (GeneInfoIT=(GeneDataIT -> second).begin() ; GeneInfoIT!=(GeneDataIT -> second).end()  ;  GeneInfoIT ++)
		{
			Start=(GeneInfoIT->second).GeneStart;
			End=(GeneInfoIT->second).GeneEnd;
			MapSSEE=RegionAAA.find(Start);
			if (MapSSEE==RegionAAA.end())
			{
				RegionAAA.insert(map <int,int>  :: value_type(Start,End)) ;
			}
			else
			{
				if (End >(MapSSEE->second))
				{
					MapSSEE->second=End;
				}
			}
		}




		MapSSEE=RegionAAA.begin() ;
		Start=MapSSEE->first;
		End=MapSSEE->second;

		map <int,int> Start2End ;
		Start2End[Start]=End;
		RegionMerger.insert(map <int,map <int,int> > ::value_type(ChrInt,Start2End));
		MergerIt=RegionMerger.find(ChrInt);
		MapSSEE++;


		for(  ; MapSSEE!=RegionAAA.end() ; MapSSEE++ )
		{
			if (  (MapSSEE->first) >  End    )
			{
				Start=MapSSEE->first;
				End=MapSSEE->second;
				(MergerIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
			}
			else if  (  (MapSSEE->second) >  End    )
			{
				MergerMapSSEE=(MergerIt->second).end(); MergerMapSSEE--;
				MergerMapSSEE->second= (MapSSEE->second) ;
				End=MapSSEE->second;
			}
		}
	}

	// check region
	if (RegionMerger.empty())
	{
		//MeMBinWindows
		int  MeMBinWindows=10000000;
		if ((paraFA04->WinSize) ==0 )
		{
			//cout<<"Warning: GFF/GTF or BED was not provided, the total chromosome length will be used as the parsing region.\n";
			(paraFA04->InInt2)=0;
		}
		else if ((paraFA04->WinSize) <150 )
		{
			(paraFA04->InInt2)=6;
		}
		else
		{
			(paraFA04->InInt2)=5;
			MeMBinWindows=(paraFA04->WinSize);
		}


		for(int i = 0; i < (header->n_targets); i++)
		{
			Start=1;
			End=2;

			for (Start=1;  End <= (header->target_len[i]) ; Start+=MeMBinWindows )
			{
				End=Start+MeMBinWindows-1;
				if (End>(header->target_len[i]))
				{
					End=(header->target_len[i]);
				}


				string GeneID=header->target_name[i]+Int2Str(Start);
				map <string,GeneInfo> TmpGene;
				GeneInfo  GeneInfoTmp ;
				GeneInfoTmp.GeneStart=Start;
				GeneInfoTmp.GeneEnd=End;
				GeneInfoTmp.GeneLength=(End-Start+1);	
				GeneInfoTmp.CDSList.push_back({Start,End});

				if (RefIn)
				{
					for (int ii=Start-1; ii<End ; ii++)
					{
						GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[i])[ii]];
					}
				}

				if (Start==1)
				{
					TmpGene[GeneID]=GeneInfoTmp;
					GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(i,TmpGene));

					map <int,int> AA;
					AA.insert(map <int,int>  :: value_type(Start,End));
					RegionMerger.insert(map <int,map <int,int> >  :: value_type(i,AA));

				}
				else
				{
					GeneDataIT=GeneData.find(i);
					(GeneDataIT->second).insert(map <string,GeneInfo>  :: value_type(GeneID,GeneInfoTmp)) ;

					MergerIt=RegionMerger.find(i);
					(MergerIt->second).insert(map <int,int>  :: value_type(Start,End)) ;
				}

				End+=2;

			}


		}

	}





	(paraFA04->InStr3)=(paraFA04->InStr3).substr(0,(paraFA04->InStr3).length()-3);
	string path=(paraFA04->InStr3);	
	ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);

	string PrefixO=path;

	if (ext == "stat" || ext == "bed")
	{
		PrefixO=path.substr(0,path.rfind('.'));
	}

	string  OutStatFile=PrefixO+".gene.stat.gz";
	string  OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";

	if  ((paraFA04->InInt2)==3)
	{
		OutStatFile=PrefixO+".bed.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tRegionID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==4)
	{
		OutStatFile=PrefixO+".bed.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==5 ||  (paraFA04->InInt2)==6 )
	{
		OutStatFile=PrefixO+".win.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==0)
	{
		OutStatFile=PrefixO+".chr.stat.gz";
		OutHeader="#Chr\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}

	if (RefIn)
	{
		RefBase.clear();
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		if  ((paraFA04->InInt2)==3)
		{
			OutHeader="#Chr\tStart\tEnd\tRegionID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==4)
		{
			OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==5|| (paraFA04->InInt2)==6)
		{
			OutHeader="#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)==0)
		{
			OutHeader="#Chr\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
	}
	ogzstream  OUT (OutStatFile.c_str());

	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<OutStatFile<<endl;
		delete  paraFA04 ; return  0;
	}
	// read index file
	string bambai=BamPath+".bai";
	string bamcsi=BamPath+".csi";
	string crambai=BamPath+".crai";
	if ( ( ( access(bambai.c_str(), 0) == 0 )  ||  (access(crambai.c_str(), 0) == 0 )   ||   (access(bamcsi.c_str(), 0) == 0 )   )  && (paraFA04->TF ) )	
	{

		if((paraFA04->SiteOutPut)  ||  ((paraFA04->InInt2)==6)  )
		{
			SiteInfo **depth = new SiteInfo *[(header->n_targets)];
			for(int i = 0; i < (header->n_targets); i++)  
			{
				int CC=(header->target_len[i])+500;
				MergerIt=RegionMerger.find(i);
				string ChrName=header->target_name[i];				
				if (MergerIt==RegionMerger.end())
				{
					CC=500;
				}
				depth[i] = new SiteInfo [CC];
				//				cerr<<ChrName<<"\t"<<CC<<endl;
				for ( int32_t j =0 ; j< CC ; j++ )
				{
					depth[i][j].Depth=0;
				}
			}





			ubit64_t GenomeLen=0;
			vector< pair<int, int> > Int2Len;
			auto RegionIt=RegionMerger.begin();
			int ChrNumRun=0;
			int MaxChrLenTen=0;
			while (RegionIt!= RegionMerger.end())
			{
				int  ChrNum=RegionIt->first;
				GenomeLen+=(header->target_len[ChrNum]);
				Int2Len.push_back({ChrNum,header->target_len[ChrNum]});
				if (header->target_len[ChrNum]>MaxChrLenTen)
				{
					MaxChrLenTen=header->target_len[ChrNum];
				}
				RegionIt++;
				ChrNumRun++;
			}


			MaxChrLenTen=int(MaxChrLenTen/10);
			RegionIt=RegionMerger.begin();
			int BigChrNum=0;
			while (RegionIt!= RegionMerger.end())
			{
				int  ChrNum=RegionIt->first;
				if (header->target_len[ChrNum]>MaxChrLenTen)
				{
					BigChrNum++;
				}
				RegionIt++;
			}

			int BamThreadSum=0;
			if ( (paraFA04->CPU) >  BigChrNum )
			{
				BamThreadSum=(paraFA04->CPU)-BigChrNum;
				(paraFA04->CPU)=BigChrNum;
			}
			ubit64_t meanLenDea=GenomeLen/(paraFA04->CPU)+2;
			map <int,vector <int> >  VecChr;
			map <int,vector <int> > :: iterator  VecChrITT;

			sort(Int2Len.begin(),Int2Len.end(),[](const pair<int, int>& a, const pair<int, int>& b){ return a.second > b.second; }); 

			int StBase=0;
			int ShiftQ=0;
			int RunCPU=(paraFA04->CPU);
			ubit64_t *CutThread = new ubit64_t [RunCPU];
			int  *BamThread =new int[RunCPU];
			for (int i = 0; i < RunCPU ;++i)
			{
				CutThread[i]=0;
				BamThread[i]=0;
			}

			while(BamThreadSum>0)
			{
				ShiftQ=StBase%RunCPU;
				StBase++;
				BamThread[ShiftQ]++;
				BamThreadSum--;
			}
			ShiftQ=0;
			StBase=0;
			for (int i = 0; i != Int2Len.size(); i++)
			{
				int  ChrNum=Int2Len[i].first;
				int  Num=int((i-ShiftQ)%RunCPU)+StBase;

				VecChrITT=VecChr.find(Num);
				if  (VecChrITT==VecChr.end())
				{
					vector <int> tmp;
					tmp.emplace_back(ChrNum);
					VecChr.insert( map <int,vector <int> >  :: value_type (Num,tmp));
				}
				else
				{
					(VecChrITT->second).emplace_back(ChrNum);
				}

				CutThread[Num]+=Int2Len[i].second;
				if (  CutThread[Num] >meanLenDea    )
				{
					StBase++;
					RunCPU--;
					ShiftQ=i+1;
				}
			}
			delete [] CutThread ;

			vector<thread> ThreadsVector ;
			VecChrITT=VecChr.begin();
			ShiftQ=0;
			while(VecChrITT!=VecChr.end())
			{
				if  (BamThread[ShiftQ]> 0 ) {BamThread[ShiftQ]=2;}
				ThreadsVector.emplace_back(ProDealChrBambaiOUTSite , std::ref(BamPath), paraFA04,std::ref(depth) ,std::ref(RegionMerger) , std::ref(VecChrITT->second),std::ref(GeneData) ,std::ref(BamThread[ShiftQ]));
				VecChrITT++;
				ShiftQ++;
			}


			int AA=ThreadsVector.size();
			for (int i=0 ; i<AA ; i++)
			{
				ThreadsVector[i].join();
			}
			ThreadsVector.clear();
			delete [] BamThread ;



			if (paraFA04->SiteOutPut)
			{
				string  OutSSiteFile=PrefixO+".SiteDepth.gz";
				ogzstream  OUTFA (OutSSiteFile.c_str());
				for(int i = 0; i < (header->n_targets); i++)
				{
					int CC=(header->target_len[i]);
					MergerIt=RegionMerger.find(i);
					if (MergerIt==RegionMerger.end())
					{
						CC=100;
						continue ;
					}
					string ChrName=header->target_name[i];
					for (int32_t j =0 ; j< CC ; j++)
					{
						OUTFA<<ChrName<<"\t"<<j<<"\t"<<depth[i][j].Depth<<"\n";
					}
				}
				OUTFA.close() ;
			}





			if (   (paraFA04->InInt2)==6 )
			{
				ubit64_t  SS_Len =0;
				ubit64_t  SS_Cov=0;
				ubit64_t  SS_TotalD=0;
				if (RefIn)
				{
					OUT<<OutHeader;
					ubit64_t SS_GCGC =0;
					for(int i = 0; i < (header->n_targets); i++)  
					{
						int CC=(header->target_len[i]);
						MergerIt=RegionMerger.find(i);
						if (MergerIt==RegionMerger.end())
						{
							CC=100;
							continue ;
						}
						string ChrName=header->target_name[i];
						int Start ; int End;
						int GeneLength ; int GeneCover ; int GeneDepth=0;
						int GeneGCGC=0;

						for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
						{

							Start=j-1; End=Start+(paraFA04->WinSize);
							if (End>CC ) {End=CC;}
							GeneCover=0;GeneDepth=0;GeneGCGC=0;

							for ( ;Start<End ; Start++)
							{
								if (depth[i][Start].Depth >= paraFA04->minDep)
								{
									GeneCover++;
									GeneDepth+=(depth[i][Start].Depth);
								}
								GeneGCGC+=GCGCArry[(RefBase[i])[Start]];
							}
							GeneLength=End-j+1;
							double Coverage=(GeneCover)*100.0/(GeneLength);
							double MeanDepth=(GeneDepth)*1.0/(GeneLength);
							double GeneGC=(GeneGCGC)*100.0/(GeneLength);
							OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\n";

							SS_Cov+=(GeneCover);
							SS_Len+=(GeneLength);
							SS_GCGC+=(GeneGCGC);
							SS_TotalD+=(GeneDepth);
						}
					}
					double Coverage=SS_Cov*100.0/SS_Len;
					double MeanDepth=SS_TotalD*1.0/SS_Len;
					double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

					OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tGC(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ALLGeneGCRation<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;


				}
				else
				{

					OUT<<OutHeader;
					for(int i = 0; i < (header->n_targets); i++)  
					{
						int CC=(header->target_len[i]);
						MergerIt=RegionMerger.find(i);
						if (MergerIt==RegionMerger.end())
						{
							CC=100;
							continue ;
						}
						string ChrName=header->target_name[i];
						int Start ; int End;
						int GeneLength ; int GeneCover ; int GeneDepth=0;

						for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
						{

							Start=j-1; End=Start+(paraFA04->WinSize);
							if (End>CC ) {End=CC;}
							GeneCover=0;GeneDepth=0;

							for ( ;Start<End ; Start++)
							{
								if (depth[i][Start].Depth >= paraFA04->minDep)
								{
									GeneCover++;
									GeneDepth+=(depth[i][Start].Depth);
								}
							}
							GeneLength=End-j+1;
							double Coverage=(GeneCover)*100.0/(GeneLength);
							double MeanDepth=(GeneDepth)*1.0/(GeneLength);
							OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<"\n";

							SS_Cov+=(GeneCover);
							SS_Len+=(GeneLength);
							SS_TotalD+=(GeneDepth);
						}
					}
					double Coverage=SS_Cov*100.0/SS_Len;
					double MeanDepth=SS_TotalD*1.0/SS_Len;

					OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;
				}
			}






			for(int i = 0; i <(header->n_targets); i++)
			{
				delete[] depth[i];  
			}
			delete[] depth;




		}
		else
		{
			ubit64_t GenomeLen=0;
			vector< pair<int, int> > Int2Len;
			auto RegionIt=RegionMerger.begin();
			int ChrNumRun=0;
			int MaxChrLenTen=0;
			while (RegionIt!= RegionMerger.end())
			{
				int  ChrNum=RegionIt->first;
				GenomeLen+=(header->target_len[ChrNum]);
				Int2Len.push_back({ChrNum,header->target_len[ChrNum]});
				if (header->target_len[ChrNum]>MaxChrLenTen)
				{
					MaxChrLenTen=header->target_len[ChrNum];
				}
				RegionIt++;
				ChrNumRun++;
			}


			MaxChrLenTen=int(MaxChrLenTen/10);
			RegionIt=RegionMerger.begin();
			int BigChrNum=0;
			while (RegionIt!= RegionMerger.end())
			{
				int  ChrNum=RegionIt->first;
				if (header->target_len[ChrNum]>MaxChrLenTen)
				{
					BigChrNum++;
				}
				RegionIt++;
			}

			int BamThreadSum=0;
			if ( (paraFA04->CPU) >  BigChrNum )
			{
				BamThreadSum=(paraFA04->CPU)-BigChrNum;
				(paraFA04->CPU)=BigChrNum;
			}
			ubit64_t meanLenDea=GenomeLen/(paraFA04->CPU)+2;
			map <int,vector <int> >  VecChr;
			map <int,vector <int> > :: iterator  VecChrITT;

			sort(Int2Len.begin(),Int2Len.end(),[](const pair<int, int>& a, const pair<int, int>& b){ return a.second > b.second; }); 

			int StBase=0;
			int ShiftQ=0;
			int RunCPU=(paraFA04->CPU);
			ubit64_t *CutThread = new ubit64_t [RunCPU];
			int  *BamThread =new int[RunCPU];
			for (int i = 0; i < RunCPU ;++i)
			{
				CutThread[i]=0;
				BamThread[i]=0;
			}

			while(BamThreadSum>0)
			{
				ShiftQ=StBase%RunCPU;
				StBase++;
				BamThread[ShiftQ]++;
				BamThreadSum--;
			}
			ShiftQ=0;
			StBase=0;
			for (int i = 0; i != Int2Len.size(); i++)
			{
				int  ChrNum=Int2Len[i].first;
				int  Num=int((i-ShiftQ)%RunCPU)+StBase;

				VecChrITT=VecChr.find(Num);
				if  (VecChrITT==VecChr.end())
				{
					vector <int> tmp;
					tmp.emplace_back(ChrNum);
					VecChr.insert( map <int,vector <int> >  :: value_type (Num,tmp));
				}
				else
				{
					(VecChrITT->second).emplace_back(ChrNum);
				}

				CutThread[Num]+=Int2Len[i].second;
				if (  CutThread[Num] >meanLenDea    )
				{
					StBase++;
					RunCPU--;
					ShiftQ=i+1;
				}
			}
			delete [] CutThread ;

			vector<thread> ThreadsVector ;
			VecChrITT=VecChr.begin();
			ShiftQ=0;
			while(VecChrITT!=VecChr.end())
			{
				if  (BamThread[ShiftQ]> 0 ) {BamThread[ShiftQ]=2;}
				ThreadsVector.emplace_back(ProDealChrBambai , std::ref(BamPath), paraFA04, std::ref(RegionMerger) , std::ref(VecChrITT->second),std::ref(GeneData) ,std::ref(BamThread[ShiftQ]));
				VecChrITT++;
				ShiftQ++;
			}

			int AA=ThreadsVector.size();
			for (int i=0 ; i<AA ; i++)
			{
				ThreadsVector[i].join();
			}
			ThreadsVector.clear();
			delete [] BamThread ;

		}




	}


	else
	{


		string headertext=header->text;
		string::size_type  posA = (headertext).find("\tSO:");
		bool sortBam=false;
		if (posA!=string::npos)
		{
			posA+=4;
			string::size_type  posB=(headertext).find_first_of("\n\t", posA);
			string  sortinfo=(headertext).substr(posA, posB-posA);
			if (sortinfo=="coordinate")
			{
				sortBam=true;
			}
		}



		SiteInfo **depth = new SiteInfo *[(header->n_targets)];
		bool  *EndChr = new bool [(header->n_targets)];
		map <int,int> :: iterator  *ArryIt = new map <int,int> :: iterator [(header->n_targets)];
		map <int,int> :: iterator  *ArryItEnd = new map <int,int> :: iterator [(header->n_targets)];


		for(int i = 0; i < (header->n_targets); i++)  
		{
			int CC=(header->target_len[i])+500;
			EndChr[i]=false ;
			MergerIt=RegionMerger.find(i);
			if (MergerIt==RegionMerger.end())
			{
				EndChr[i]=true ;
				if (sortBam) { CC=500;}
			}
			else
			{
				ArryIt[i]=(MergerIt->second).begin();
				ArryItEnd[i]=(MergerIt->second).end();
			}

			depth[i] = new SiteInfo  [CC];

			for (int32_t j =0 ; j< CC ; j++)
			{
				depth[i][j].Depth=0;
			}
		}

		bam1_t *aln = bam_init1();
		uint32_t *cigar;

		uint32_t flags = (paraFA04->flags);

		bam_hdr_t *headerRR;
		samFile *BamInRR = hts_open(BamPath.c_str(), "r");
		hts_set_log_level(HTS_LOG_OFF);
		int numThreads = (paraFA04->CPU);
		hts_set_opt(BamInRR, HTS_OPT_NTHREADS, numThreads);
		if(BamInRR->format.format == htsExactFormat::cram)
		{	
			const char* ref_file = (paraFA04->reference).c_str();
			hts_set_fai_filename(BamInRR, ref_file);
			hts_set_opt(BamInRR, CRAM_OPT_DECODE_MD, 0);
			hts_set_opt(BamInRR, CRAM_OPT_REQUIRED_FIELDS, SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR);

		}

		headerRR = sam_hdr_read(BamInRR);

		if (sortBam)
		{
			cout <<"Warning: PanDepth will run in No Index mode: "<<BamPath<<endl;
			bool ALLEnd=false;
			while (sam_read1(BamInRR, header, aln) >= 0)
			{
				if ( EndChr[(aln->core).tid])  {continue ;}	
				if ( (aln->core).qual < (paraFA04->InInt) )
				{
					continue ;
				}
				if ( aln->core.flag & flags ) continue;
				int32_t EndRead=bam_endpos(aln);
				if (EndRead <(ArryIt[(aln->core).tid]->first)) {continue;}

				int32_t StartRead=((aln->core).pos);
				if ( StartRead  > (ArryIt[(aln->core).tid]->second) )
				{

					for (ArryIt[(aln->core).tid]++  ; ArryIt[(aln->core).tid]!=ArryItEnd[(aln->core).tid] ;ArryIt[(aln->core).tid]++)
					{
						if ( StartRead <=  (ArryIt[(aln->core).tid]->second)  )
						{
							break;
						}
					}
					if  (ArryIt[(aln->core).tid]==ArryItEnd[(aln->core).tid])
					{
						EndChr[(aln->core).tid]=true;
						ALLEnd=true;
						for(int i = 0; i < (header->n_targets); i++)
						{
							if (!EndChr[i])
							{
								ALLEnd=false;
							}
						}
						if (ALLEnd)
						{
							break;
						}
					}
				}

				cigar = bam_get_cigar(aln);
				int  endTmp ;
				for(int i=0; i < aln->core.n_cigar;++i)
				{
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					switch (cig)
					{
						case 0:
						case 7:
						case 8 :
							endTmp=StartRead+ncig;
							for (  ; StartRead<endTmp;StartRead++)
							{
								(depth[(aln->core).tid][StartRead]).Depth++;
							}
							break;
						case 2:
						case 3:
							StartRead=StartRead+ncig;
							break;
					}
				}
			}


		}
		else
		{
			cout <<"Warning: Can't find index file of input BAM/CRAM. PanDepth will run in No Index mode: "<<BamPath<<endl;
			bool ALLEnd=false;
			while (sam_read1(BamInRR, header, aln) >= 0)
			{
				if ( (aln->core).qual < (paraFA04->InInt) )
				{
					continue ;
				}
				if ( aln->core.flag & flags ) continue;
				int32_t EndRead=bam_endpos(aln);
				int32_t StartRead=((aln->core).pos);
				cigar = bam_get_cigar(aln);
				int  endTmp ;
				for(int i=0; i < aln->core.n_cigar;++i)
				{
					int cig=bam_cigar_op(cigar[i]);
					int ncig = bam_cigar_oplen(cigar[i]);
					switch (cig)
					{
						case 0:
						case 7:
						case 8 :
							endTmp=StartRead+ncig;
							for (  ; StartRead<endTmp;StartRead++)
							{
								(depth[(aln->core).tid][StartRead]).Depth++;
							}
							break;
						case 2:
						case 3:
							StartRead=StartRead+ncig;
							break;
					}
				}
			}

		}


		sam_close(BamInRR);
		bam_hdr_destroy(headerRR);

		for(int i = 0; i < (header->n_targets); i++)  
		{
			MergerIt=RegionMerger.find(i);
			if (MergerIt==RegionMerger.end())
			{

			}
			else
			{
				GeneDataIT=GeneData.find(i);
				StatChrDepthLowMEM ( depth[i] , GeneDataIT->second, paraFA04->minDep);
			}
		}

		bam_destroy1(aln);


		if(paraFA04->SiteOutPut)
		{
			string  OutSSiteFile=PrefixO+".SiteDepth.gz";
			ogzstream  OUTFA (OutSSiteFile.c_str());
			for(int i = 0; i < (header->n_targets); i++)  
			{
				int CC=(header->target_len[i]);
				MergerIt=RegionMerger.find(i);
				if (MergerIt==RegionMerger.end())
				{
					CC=100;
					continue ;
				}
				string ChrName=header->target_name[i];
				for (int32_t j =0 ; j< CC ; j++)
				{
					OUTFA<<ChrName<<"\t"<<j<<"\t"<<depth[i][j].Depth<<"\n";
				}
			}
			OUTFA.close() ;
		}


		if (   (paraFA04->InInt2)==6 )
		{
			ubit64_t  SS_Len =0;
			ubit64_t  SS_Cov=0;
			ubit64_t  SS_TotalD=0;
			if (RefIn)
			{
				OUT<<OutHeader;
				ubit64_t SS_GCGC =0;
				for(int i = 0; i < (header->n_targets); i++)  
				{
					int CC=(header->target_len[i]);
					MergerIt=RegionMerger.find(i);
					if (MergerIt==RegionMerger.end())
					{
						CC=100;
						continue ;
					}
					string ChrName=header->target_name[i];
					int Start ; int End;
					int GeneLength ; int GeneCover ; int GeneDepth=0;
					int GeneGCGC=0;

					for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
					{

						Start=j-1; End=Start+(paraFA04->WinSize);
						if (End>CC ) {End=CC;}
						GeneCover=0;GeneDepth=0;GeneGCGC=0;

						for ( ;Start<End ; Start++)
						{
							if (depth[i][Start].Depth >= paraFA04->minDep)
							{
								GeneCover++;
								GeneDepth+=(depth[i][Start].Depth);
							}
							GeneGCGC+=GCGCArry[(RefBase[i])[Start]];
						}
						GeneLength=End-j+1;
						double Coverage=(GeneCover)*100.0/(GeneLength);
						double MeanDepth=(GeneDepth)*1.0/(GeneLength);
						double GeneGC=(GeneGCGC)*100.0/(GeneLength);
						OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\n";

						SS_Cov+=(GeneCover);
						SS_Len+=(GeneLength);
						SS_GCGC+=(GeneGCGC);
						SS_TotalD+=(GeneDepth);
					}
				}
				double Coverage=SS_Cov*100.0/SS_Len;
				double MeanDepth=SS_TotalD*1.0/SS_Len;
				double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

				OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tGC(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ALLGeneGCRation<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;


			}
			else
			{

				OUT<<OutHeader;
				for(int i = 0; i < (header->n_targets); i++)  
				{
					int CC=(header->target_len[i]);
					MergerIt=RegionMerger.find(i);
					if (MergerIt==RegionMerger.end())
					{
						CC=100;
						continue ;
					}
					string ChrName=header->target_name[i];
					int Start ; int End;
					int GeneLength ; int GeneCover ; int GeneDepth=0;

					for (int32_t j =1 ; j< CC ; j+=((paraFA04->WinSize)))
					{

						Start=j-1; End=Start+(paraFA04->WinSize);
						if (End>CC ) {End=CC;}
						GeneCover=0;GeneDepth=0;

						for ( ;Start<End ; Start++)
						{
							if (depth[i][Start].Depth >= paraFA04->minDep)
							{
								GeneCover++;
								GeneDepth+=(depth[i][Start].Depth);
							}
						}
						GeneLength=End-j+1;
						double Coverage=(GeneCover)*100.0/(GeneLength);
						double MeanDepth=(GeneDepth)*1.0/(GeneLength);
						OUT<<ChrName<<"\t"<<j<<"\t"<<End<<"\t"<<GeneLength<<"\t"<<GeneCover<<"\t"<<(GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<"\n";

						SS_Cov+=(GeneCover);
						SS_Len+=(GeneLength);
						SS_TotalD+=(GeneDepth);
					}
				}
				double Coverage=SS_Cov*100.0/SS_Len;
				double MeanDepth=SS_TotalD*1.0/SS_Len;

				OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;
			}
		}




		for(int i = 0; i <(header->n_targets); i++)
		{
			delete[] depth[i];  
		}
		delete[] depth;
		delete[] EndChr;
		delete[] ArryItEnd;
	}

	cout <<"INFO: Input data read done"<<endl;

	ubit64_t  SS_Len =0;
	ubit64_t  SS_Cov=0;
	ubit64_t  SS_TotalD=0;

	if ((RefIn)  &&  ((paraFA04->InInt2)!=6))
	{
		OUT<<OutHeader;
		ubit64_t  SS_GCGC=0;
		if  ((paraFA04->InInt2)==1 ||  (paraFA04->InInt2)==2   ||  (paraFA04->InInt2)==3   ||  (paraFA04->InInt2)==4 )
		{
			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 
				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);
					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<GeneID<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}
				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		else if ((paraFA04->InInt2)==0)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];

				ubit64_t  ChrSS_Len =0;
				ubit64_t  ChrSS_Cov=0;
				ubit64_t  ChrSS_TotalD=0;
				ubit64_t  ChrSS_GCGC=0;

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					ChrSS_Cov+=((GeneInfoIT->second).GeneCover);
					ChrSS_Len+=((GeneInfoIT->second).GeneLength);
					ChrSS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					ChrSS_TotalD+=((GeneInfoIT->second).GeneDepth);
				}

				double GeneGC=ChrSS_GCGC*100.0/ChrSS_Len;
				double MeanDepth=ChrSS_TotalD*1.0/ChrSS_Len;
				double Coverage=ChrSS_Cov*100.0/ChrSS_Len;
				OUT<<ChrName<<"\t"<<ChrSS_Len<<"\t"<<ChrSS_Cov<<"\t"<<ChrSS_TotalD<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth<<"\n";

			}
		}
		else  if ((paraFA04->InInt2)==5)
		{


			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_GCGC+=((GeneInfoIT->second).GeneGCGC);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		double Coverage=SS_Cov*100.0/SS_Len;
		double MeanDepth=SS_TotalD*1.0/SS_Len;
		double ALLGeneGCRation=SS_GCGC*100.0/SS_Len;

		OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tGC(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<ALLGeneGCRation<<"\tCoverage(%): "<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;

	}
	else if (((paraFA04->InInt2)!=6))
	{
		OUT<<OutHeader;
		if  ((paraFA04->InInt2)==1 ||  (paraFA04->InInt2)==2   ||  (paraFA04->InInt2)==3   ||  (paraFA04->InInt2)==4 )
		{
			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					double GeneGC=((GeneInfoIT->second).GeneGCGC)*100.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);
					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<GeneID<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth;

					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}

				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		else if ((paraFA04->InInt2)==0)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];

				ubit64_t  ChrSS_Len =0;
				ubit64_t  ChrSS_Cov=0;
				ubit64_t  ChrSS_TotalD=0;

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					ChrSS_Cov+=((GeneInfoIT->second).GeneCover);
					ChrSS_Len+=((GeneInfoIT->second).GeneLength);
					ChrSS_TotalD+=((GeneInfoIT->second).GeneDepth);
				}

				double MeanDepth=ChrSS_TotalD*1.0/ChrSS_Len;
				double Coverage=ChrSS_Cov*100.0/ChrSS_Len;
				OUT<<ChrName<<"\t"<<ChrSS_Len<<"\t"<<ChrSS_Cov<<"\t"<<ChrSS_TotalD<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth<<"\n";

			}
		}
		else if ((paraFA04->InInt2)==5)
		{

			for ( GeneDataIT=GeneData.begin();  GeneDataIT != GeneData.end(); ++GeneDataIT)
			{
				int  ChrNum=GeneDataIT->first ;
				string ChrName=header->target_name[ChrNum];
				map <int , string >   SortResult ;
				map <int, string  > :: iterator   SortResultIT; 

				for (GeneInfoIT=(GeneDataIT->second).begin() ; GeneInfoIT!=(GeneDataIT->second).end(); GeneInfoIT++)
				{
					string  GeneID=GeneInfoIT->first;
					double Coverage=((GeneInfoIT->second).GeneCover)*100.0/((GeneInfoIT->second).GeneLength);
					double MeanDepth=((GeneInfoIT->second).GeneDepth)*1.0/((GeneInfoIT->second).GeneLength);
					SS_Cov+=((GeneInfoIT->second).GeneCover);
					SS_Len+=((GeneInfoIT->second).GeneLength);
					SS_TotalD+=((GeneInfoIT->second).GeneDepth);

					stringstream  ss ;
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneStart)<<"\t"<<((GeneInfoIT->second).GeneEnd)<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth;
					SortResultIT=SortResult.find(((GeneInfoIT->second).GeneStart));
					if (SortResultIT==SortResult.end())
					{
						SortResult.insert( map <int,string>  :: value_type(((GeneInfoIT->second).GeneStart),ss.str()));
					}
					else
					{
						(SortResultIT->second)+="\n";
						(SortResultIT->second)+=ss.str();
					}
				}

				for(SortResultIT=SortResult.begin(); SortResultIT!=SortResult.end(); SortResultIT++)
				{
					OUT<<(SortResultIT->second) <<"\n";
				}
			}
		}
		double Coverage=SS_Cov*100.0/SS_Len;
		double MeanDepth=SS_TotalD*1.0/SS_Len;

		OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tCoverage(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;

	}


	OUT.close();


	bam_hdr_destroy(header);
	delete paraFA04 ;

	return 0;
}



#endif






