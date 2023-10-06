/* The MIT License

Copyright (c) 2023- by Huiyang Yu,
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


using namespace std;
typedef unsigned long long ubit64_t;
//KSEQ_INIT(gzFile, gzread)


void bamCov_help()
{
	cout<<""
		"Usage: pandepth -i in.bam [-g gene.gff|-b region.bed] -o  outPrefix\n"
		" Input/Output options\n"
		"   -i    <str>     input of bam/cram file\n"
		"   -o    <str>     prefix of output file\n"
		" Target options\n"
		"   -g    <str>     input gff/gtf file for gene region\n"
		"   -f    <str>     gff/gtf feature type to parse, CDS or exon [CDS]\n"
		"   -b    <str>     input bed file for list of regions\n"
		" Filter options\n"
		"   -q    <int>     min mapping quality [0]\n"
		"   -x    <int>     exclude reads with any of the bits in FLAG set [1796]\n"
		" Other options\n"
		"   -t    <int>     number of threads [3]\n"
//		"   -l    <int>     the region extension length [auto]\n"
		"   -r    <str>     input reference genome file for cram decode or GC parse\n"
		"   -c              enable the calculation of GC content within the target region (requires -r)\n"
		"   -h              show this help [v2.16]\n"
		"\n";
}

int bamCov_help01(int argc, char **argv , In3str1v * paraFA04   )
{
	if (argc <=1 ) {bamCov_help();return 0;}
	int file_count=0;

	//	bool AA=true;	bool BB=true;	bool CC=true;
	//paraFA04->flags = 1796 ;	
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
			paraFA04->InStr1=argv[i];
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
		else if (flag == "l")
		{
			if(i + 1 == argc) {LogLackArg(flag); return 0;}
			i++;
			paraFA04->readoverLen=atoi(argv[i]);
			if ( ( paraFA04->readoverLen) >=0)  {paraFA04->readoverLen=0;}
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
	//

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
	return 1 ;
}
//
void  StatChrDepth ( short unsigned int *depth  ,  map <string,GeneInfo>  &  GeneStat )
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
				if ( depth[ii]>0 )
				{
					(iter->second).GeneCover++;
					(iter->second).GeneDepth+=depth[ii];
				}
			}
		}

	}
}
//


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
		for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++ )
		{
			long beg=(MapSSEE->first);
			if  (beg<1) {beg=1;}
			long end=(MapSSEE->second);
			if (end>ChrLen) {end=ChrLen;}
			RegionArry[CountRegion++] = CharMap;
			CharMap += (sprintf(CharMap, "%s:%lu-%lu", ChrName.c_str(), beg, end)+1);
		}

		ChrLen+=500;
		unsigned short int *depth = new unsigned short int [ChrLen];
		for (int  ccv=0 ; ccv<ChrLen  ; ccv++)
		{
			depth[ccv]=0;
		}





		bam_mplp_t mplp;
		int	pos;

		rcnt = CountRegion;
		iter = sam_itr_regarray(idx, headerAA, RegionArry, rcnt);
//		cerr<<RegionArry[0]<<"\t"<<NumberRegionSize<<"\t"<<rcnt<<endl;
	

		delete [] RegionArry ;
//		delete CharMap ;

		//while ((ret = sam_itr_multi_next(fphts, iter, aln)) >= 0)
		while ((ret = sam_itr_next(fphts, iter, aln)) >= 0)
		{
			if ( aln->core.flag & flags )  { continue; }
			if ( (aln->core).qual < (paraFA04->InInt) )	{	continue ;	}

			cigar = bam_get_cigar(aln);
			int32_t StartRead=((aln->core).pos);
			int32_t endTmp;
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
		StatChrDepth ( depth  , GeneDataIT->second );

		delete [] depth;

	}

	bam_destroy1(aln);
	if (idx) hts_idx_destroy(idx);
	if (iter) sam_itr_destroy(iter);
	if (headerAA) sam_hdr_destroy(headerAA);
}

//int bamCov_main(int argc, char *argv[])
int main(int argc, char *argv[])
{
	In3str1v *paraFA04 = new In3str1v;
	paraFA04->InInt=-1;
	//paraFA04->TF=false ;
	if ((bamCov_help01(argc, argv, paraFA04)==0))
	{
		delete paraFA04 ;
		return 0 ;
	}

	string  BamPath=(paraFA04->InStr1);
	if (BamPath.length()<=0)  
	{
		cerr<<"Error: Failed to open the BAM/CRAM file: "<<BamPath<<endl;
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
	int numThreads = 1;
	hts_set_opt(BamIn, HTS_OPT_NTHREADS, numThreads);
	header = sam_hdr_read(BamIn);

	map <string,int> Chr2IntMap; 
	for(int i = 0; i < (header->n_targets); i++) 
	{
		string ChrName=header->target_name[i];
		Chr2IntMap.insert( map <string,int>  :: value_type (ChrName,i));
	}

	bam1_t *alnTA = bam_init1();
	int readoverLen=0;
/*
	int FlagReadCount=0;
	while ( (sam_read1(BamIn, header, alnTA) >= 0)  && (FlagReadCount<1688) )
	{
		FlagReadCount++;
		if ( readoverLen < ((alnTA->core).l_qseq))
		{
			readoverLen=((alnTA->core).l_qseq);
		}
	}
	bam_destroy1(alnTA);
	if (readoverLen>3000  && readoverLen<100000 ) { readoverLen=100000;}
	if ((paraFA04->readoverLen)>=150)
	{
		readoverLen=paraFA04->readoverLen;
	}
*/
	sam_close(BamIn);

	///////// swimming in the sky and flying in the sea ////////////
	//	/*
	map <int,string>  RefBase ;
	bool  RefIn=false;
	if ((paraFA04->gc)==true){
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
	//	cerr<<"read ref.fa done"<<endl;

	map <int,map <int,int> > Region;
	map <int,map <int,int> > :: iterator  RegionIt ;
	map <int,int> :: iterator MapSSEE ;

	map <string,int> :: iterator  MapItChr2Int ;
	unordered_map <string,int> :: iterator  UnMapItChr2Int ;
	int Start ; int End ;  string ChrName ;
	string GeneID ;
	map <int,map <string,GeneInfo> > GeneData;
	map <int,map <string,GeneInfo> >  :: iterator  GeneDataIT;
	map <string,GeneInfo>   :: iterator  GeneInfoIT  ;

	//gff  1    /		gtf  2  exon  /  bed 3  /  bed 4   / examp  0

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
						Start--;
						if (RefIn)
						{
							for (int ii=Start; ii<End ; ii++)
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
							Start--;
							if (RefIn)
							{
								for (int ii=Start; ii<End ; ii++)
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

					Start=Start-readoverLen;End=End+readoverLen;
					RegionIt=Region.find(MapItChr2Int->second);
					if (RegionIt==Region.end())
					{
						map <int,int> Start2End ;
						Start2End[Start]=End;
						Region.insert(map <int,map <int,int> > ::value_type(MapItChr2Int->second,Start2End));
					}
					else
					{
						MapSSEE=(RegionIt->second).find(Start);
						if (MapSSEE==(RegionIt->second).end())
						{
							(RegionIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
						}
						else
						{
							if  (End>(MapSSEE->second))
							{
								MapSSEE->second=End;
							}
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
						Start--;
						if (RefIn)
						{
							for (int ii=Start; ii<End ; ii++)
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
							Start--;
							if (RefIn)
							{
								for (int ii=Start; ii<End ; ii++)
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

					Start=Start-readoverLen;End=End+readoverLen;
					RegionIt=Region.find(MapItChr2Int->second);
					if (RegionIt==Region.end())
					{
						map <int,int> Start2End ;
						Start2End[Start]=End;
						Region.insert(map <int,map <int,int> > ::value_type(MapItChr2Int->second,Start2End));
					}
					else
					{
						MapSSEE=(RegionIt->second).find(Start);
						if (MapSSEE==(RegionIt->second).end())
						{
							(RegionIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
						}
						else
						{
							if  (End>(MapSSEE->second))
							{
								MapSSEE->second=End;
							}
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
						Start--;
						if (RefIn)
						{
							for (int ii=Start; ii<End ; ii++)
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
							Start--;
							if (RefIn)
							{
								for (int ii=Start; ii<End ; ii++)
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

					Start=Start-readoverLen;End=End+readoverLen;
					RegionIt=Region.find(MapItChr2Int->second);
					if (RegionIt==Region.end())
					{
						map <int,int> Start2End ;
						Start2End[Start]=End;
						Region.insert(map <int,map <int,int> > ::value_type(MapItChr2Int->second,Start2End));
					}
					else
					{
						MapSSEE=(RegionIt->second).find(Start);
						if (MapSSEE==(RegionIt->second).end())
						{
							(RegionIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
						}
						else
						{
							if  (End>(MapSSEE->second))
							{
								MapSSEE->second=End;
							}
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
						Start--;
						if (RefIn)
						{
							for (int ii=Start; ii<End ; ii++)
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
							Start--;
							if (RefIn)
							{
								for (int ii=Start; ii<End ; ii++)
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

					Start=Start-readoverLen;End=End+readoverLen;
					RegionIt=Region.find(MapItChr2Int->second);
					if (RegionIt==Region.end())
					{
						map <int,int> Start2End ;
						Start2End[Start]=End;
						Region.insert(map <int,map <int,int> > ::value_type(MapItChr2Int->second,Start2End));
					}
					else
					{
						MapSSEE=(RegionIt->second).find(Start);
						if (MapSSEE==(RegionIt->second).end())
						{
							(RegionIt->second).insert(map <int,int>  :: value_type(Start,End)) ;	
						}
						else
						{
							if  (End>(MapSSEE->second))
							{
								MapSSEE->second=End;
							}
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

	RegionIt=Region.begin() ;  	MapSSEE=(RegionIt->second).begin() ;

	for(  RegionIt=Region.begin() ; RegionIt!= Region.end(); RegionIt++ )
	{
		MapSSEE=(RegionIt->second).begin() ;
		int ChrInt = RegionIt->first ;
		Start=MapSSEE->first;
		End=MapSSEE->second;

		map <int,int> Start2End ;
		Start2End[Start]=End;
		RegionMerger.insert(map <int,map <int,int> > ::value_type(ChrInt,Start2End));
		MergerIt=RegionMerger.find(ChrInt);

		MapSSEE++;


		for(  ; MapSSEE!=(RegionIt->second).end() ; MapSSEE++ )		
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
		for(int i = 0; i < (header->n_targets); i++)  
		{
			map <int,int> AA;
			Start=1;
			End=(header->target_len[i]);
			AA.insert(map <int,int>  :: value_type(Start,End));
			RegionMerger.insert(map <int,map <int,int> >  :: value_type(i,AA));

			string GeneID=header->target_name[i];

			map <string,GeneInfo> TmpGene;
			GeneInfo  GeneInfoTmp ;
			GeneInfoTmp.GeneStart=1;
			GeneInfoTmp.GeneEnd=End;
			GeneInfoTmp.GeneLength=End;	
			GeneInfoTmp.CDSList.push_back({1,End});

			if (RefIn)
			{
				for (int ii=0; ii<End ; ii++)
				{
					GeneInfoTmp.GeneGCGC+=GCGCArry[(RefBase[i])[ii]];
				}
			}
			TmpGene[GeneID]=GeneInfoTmp;
			GeneData.insert( map <int,map <string,GeneInfo> > ::  value_type(i,TmpGene));

		}
		cout<<"Warning: GFF/GTF or BED was not provided, the total chromosome length will be used as the parsing region.\n";
		(paraFA04->InInt2)=0;

	}
	else
	{
		Region.clear();
	}


	if (RefIn)
	{
		RefBase.clear();
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
	else if ((paraFA04->InInt2)>=4)
	{
		OutStatFile=PrefixO+".bed.stat.gz";
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}
	else if ((paraFA04->InInt2)==0)
	{
		OutStatFile=PrefixO+".chr.stat.gz";
		OutHeader="#Chr\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n";
	}

	if (RefIn  &&  (paraFA04->gc) )
	{
		OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";

		if  ((paraFA04->InInt2)==3)
		{
			OutHeader="#Chr\tStart\tEnd\tRegionID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
		}
		else if ((paraFA04->InInt2)>=4)
		{
			OutHeader="#Chr\tStart\tEnd\tGeneID\tLength\tCoveredSite\tTotalDepth\tGC(%)\tCoverage(%)\tMeanDepth\n";
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
	string crambai=BamPath+".crai";
	if ( ( ( access(bambai.c_str(), 0) == 0 )  ||  (access(crambai.c_str(), 0) == 0 )  )  && (paraFA04->TF ) )	
	{

		ubit64_t GenomeLen=0;
		vector< pair<int, int> > Int2Len;
		RegionIt=RegionMerger.begin();
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


	else
	{
		unsigned short int **depth = new unsigned short int*[(header->n_targets)];
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
				CC=500;
			}
			else
			{
				ArryIt[i]=(MergerIt->second).begin();
				ArryItEnd[i]=(MergerIt->second).end();
			}

			depth[i] = new unsigned short int [CC];

			for (int32_t j =0 ; j< CC ; j++)
			{
				depth[i][j]=0;
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
		//
		headerRR = sam_hdr_read(BamInRR);

		cout <<"Warning: Cannot find index file of input BAM/CRAM. PanDepth will run in No Index mode: "<<BamPath<<endl;
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
							depth[(aln->core).tid][StartRead]++;
						}
						break;
					case 2:
					case 3:
						StartRead=StartRead+ncig;
						break;
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
				StatChrDepth ( depth[i] , GeneDataIT->second );
			}
		}

		bam_destroy1(aln);

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

	if (RefIn)
	{
		OUT<<OutHeader;
		ubit64_t  SS_GCGC=0;
		if  ((paraFA04->InInt2)!=0)
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
		else
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
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<GeneGC<<"\t"<<Coverage<<"\t"<<MeanDepth;
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
	else
	{
		OUT<<OutHeader;
		if  ((paraFA04->InInt2)!=0)
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
		else
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
					ss<<ChrName<<"\t"<<((GeneInfoIT->second).GeneLength)<<"\t"<<((GeneInfoIT->second).GeneCover)<<"\t"<<((GeneInfoIT->second).GeneDepth)<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\t"<<MeanDepth;
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

		double Coverage=SS_Cov*100.0/SS_Len ;
		double MeanDepth=SS_TotalD*1.0/SS_Len;

		OUT<<"##RegionLength: "<<SS_Len<<"\tCoveredSite: "<<SS_Cov<<"\tCoverage(%): "<<setiosflags(ios::fixed)<<setiosflags(ios::right)<<setprecision(2)<<Coverage<<"\tMeanDepth: "<<MeanDepth<<endl;

	}


	OUT.close();


	bam_hdr_destroy(header);
	delete paraFA04 ;

	return 0;
}
#endif





