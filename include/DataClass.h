
/*
 * DataClass.h
 */

#ifndef DataClass_H_
#define DataClass_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "comm.h"

using namespace std;
typedef unsigned long long ubit64_t;


///////////////// q_seq  for site ////////


class In3str1v {
	public:
		bool TF ;
		int  InInt ;
		int  InInt2 ;
		int  CPU ;
		bool TF2 ;
		bool gc ;
		bool SiteOutPut ;
		uint32_t flags ;
		int WinSize;
		int minDep;
		string InStr1 ;
		vector <string> InStr1List ;
		string InStr2 ;
		string InStr3 ;
		string CDS ;
		string reference ;
		vector <string> List ;
		In3str1v()
		{
			InStr1="";
			InStr2="";
			InStr3="";
			CDS="CDS";
			reference="";
			gc=false;
			TF=true ;
			TF2=true ;
			SiteOutPut=false;
			InInt=0 ;
			InInt2=0;
			flags=1796;
			CPU=3;
			WinSize=0;
			minDep=1;
		}
};


class  GeneInfo {
	public:
		int GeneCover;
		ubit64_t GeneDepth;
		int GeneStart;
		int GeneEnd;
		int GeneGCGC;
		ubit64_t GeneLength=0;
		vector< pair<int, int> > CDSList; 
		GeneInfo()
		{
			GeneCover=0;
			GeneDepth=0;
			GeneStart=0;
			GeneEnd=0;
			GeneLength=0;
			GeneGCGC=0;
		}
};



struct SiteInfo
{
	   unsigned int Depth: 18;
};



#endif /* DataClass_H_ */

