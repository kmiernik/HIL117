#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include <bitset>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <map>
#include <vector>
#include <TArrayS.h>



#include <stdio.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


template<std::size_t R, std::size_t L, std::size_t N>
std::bitset<N> project_range(std::bitset<N> b)
{
  static_assert(R <= L && L <= N, "invalid bitrange");
  b = b << ( N - L - 1);            // drop R rightmost bits
  b = b >> (R + (N - L - 1));  // drop L-1 leftmost bits
  return b;
}

#define dbug_line 	std::cout << "I am at line " <<  __LINE__ << std::endl;
bool debug = false;
struct dataFormatPSD{
  std::bitset<32> dataFormat;
  void     SetDataFormat(uint32_t form){dataFormat = form;}
  bool     enaDT()         const {return dataFormat.test(31);}
  bool     enaCharge()     const {return dataFormat.test(30);}
  bool     enaTS()         const {return dataFormat.test(29);}
  bool     enaExtras()     const {return dataFormat.test(28);}
  bool     enaTrace()      const {return dataFormat.test(27);}
  uint8_t  confExtras()    const {return static_cast<uint8_t>(project_range<24,26,32>(dataFormat).to_ulong());}
  uint8_t  confAnaProbe()  const {return static_cast<uint8_t>(project_range<22,23,32>(dataFormat).to_ulong());}
  uint8_t  confDigProbe1() const {return static_cast<uint8_t>(project_range<19,21,32>(dataFormat).to_ulong());}
  uint8_t  confDigProbe2() const {return static_cast<uint8_t>(project_range<16,18,32>(dataFormat).to_ulong());}
  uint16_t numSamples()    const {return static_cast<uint16_t>(project_range<0,15,32>(dataFormat).to_ulong())*8;}
  uint8_t evtSize() const{
    uint8_t evtSize = 0;
    //    if(enaDT())     evtSize+=1;
    if(enaCharge()) evtSize+=1;
    if(enaTS())     evtSize+=1;
    if(enaExtras()) evtSize+=1;
    if(enaTrace())  evtSize += numSamples()/2;
    return evtSize;
  } 
  void printDataFormat(){
    
    printf("Data format (%08X): ",dataFormat.to_ulong());
    std::cout << dataFormat << std::endl;
    printf(" --- Dual Trace [31] are %s\n",((enaDT())? "enabled" : "disabled"));
    printf(" --- Charge [30] are %s\n",((enaCharge())? "enabled" : "disabled"));
    printf(" --- TimeTag [29] are %s\n",((enaTS())? "enabled" : "disabled"));
    printf(" --- Extras [28] are %s\n",((enaExtras())? "enabled" : "disabled"));
    if(enaExtras())
      printf(" ---- Extras configuration [26:24] set to %i \n",confExtras()); 
    printf(" --- Trace [27] are %s\n",((enaTrace())? "enabled" : "disabled"));       
    printf(" --- Event length %i 32-bit words\n",evtSize());
  }
} dataForm;

struct agataKey {
  // POD with 5 uint32_t
  uint32_t size  = 0;  // in bytes
  uint32_t key   = 0;   //
  uint32_t evnum = 0;   // AGATA: defined but never used === GALILEO: domain_number/crystal_id info for \
  the EventBuilder
  uint32_t ts_0 = 0;    //
  uint32_t ts_1 = 0;    //

  uint32_t GetBytes () const { return size; }
  uint32_t GetEvnum () const { return  evnum; }
  uint64_t GetTstamp() const { return (uint64_t(ts_1 )<<32) | (uint64_t(ts_0));}
  //agataKey() { memset(this, 0, sizeof(agataKey)); }
};

struct subDataPHA_t
{
  uint16_t energy = 0 ;
  uint16_t cfd    = 0 ;   // in ps
};    

struct subDataPSD_t
{
  uint16_t qshort = 0.;
  uint16_t qlong  = 0.;
  float cfd       = 0.;   // in 10 ns units
};
struct caenEventPSD_t
{
  agataKey     theKey ;
  subDataPSD_t theData;
};

struct caenEventPHA_t
{
  agataKey     theKey ;
  subDataPHA_t theData;
};
    
caenEventPSD_t *outEvt;

struct dataFormatPHA{
  std::bitset<32> dataFormat;
  void  SetDataFormat(uint32_t form){dataFormat = form;}
  bool  enaDT()      const {return dataFormat.test(31);}
  bool  enaCharge()  const {return dataFormat.test(30);}
  bool  enaTS()      const {return dataFormat.test(29);}
  bool  enaExtras()  const {return dataFormat.test(28);}
  bool  enaTrace()   const {return dataFormat.test(27);}
  uint8_t  confExtras()    const {return static_cast<uint8_t>(project_range<24,26,32>(dataFormat).to_ulong());}
  uint8_t  confAnaProbe()  const {return static_cast<uint8_t>(project_range<22,23,32>(dataFormat).to_ulong());}
  uint8_t  confDigProbe1() const {return static_cast<uint8_t>(project_range<19,21,32>(dataFormat).to_ulong());}
  uint8_t  confDigProbe2() const {return static_cast<uint8_t>(project_range<16,18,32>(dataFormat).to_ulong());}
  uint16_t numSamples()    const {return static_cast<uint16_t>(project_range<0,15,32>(dataFormat).to_ulong())*8;}
  uint8_t evtSize() const{
    uint8_t evtSize = 0;
    if(enaCharge()) evtSize+=1;
    if(enaTS())     evtSize+=1;
    if(enaExtras()) evtSize+=1;
    if(enaTrace())  evtSize += numSamples()/2;
    return evtSize;
  }
  void printDataFormat(){

    printf("Data format (%08X): ",dataFormat.to_ulong());
    std::cout << dataFormat << std::endl;
    printf(" --- Dual Trace [31] are %s\n",((enaDT())? "enabled" : "disabled"));
    printf(" --- Charge [30] are %s\n",((enaCharge())? "enabled" : "disabled"));
    printf(" --- TimeTag [29] are %s\n",((enaTS())? "enabled" : "disabled"));
    printf(" --- Extras [28] are %s\n",((enaExtras())? "enabled" : "disabled"));
    if(enaExtras())
      printf(" ---- Extras configuration [26:24] set to %i \n",confExtras());
    printf(" --- Trace [27] are %s\n",((enaTrace())? "enabled" : "disabled"));
    printf(" --- Event length %i 32-bit words\n",evtSize());
  }
} dataFormPHA;


          


int main(int argc, char* argv[])
{
  if(argc < 3){
    std::cerr << "Usage: program InputFileName OutputFileName [maxAggr]" << std::endl;
    return -1;
  }
  long nbMaxAggr = 0;
  if(argc==4)nbMaxAggr = std::atol(argv[3]);
  FILE* fInput = fopen(argv[1],"rb");
  if(fInput == nullptr){
    std::cerr << "Error! Cannot open input file, please check that " << argv[1]
	      << " exists!" << std::endl;
    return -2;
  }
  
  TFile *fOutput = new TFile(argv[2],"recreate");
  TTree *fTree   = new TTree("rawCaen","rawCaen");


  uint16_t boardId(0), channelId(0);
  uint16_t energy(0), time(0);
  uint16_t qshort(0), qlong(0);
  uint64_t tstamp, tstamp_10ns;
  bool tsat(0), isat(0), pur(0);
  float cfd=0.;
  bool TreeEna = true;

  std::vector<uint16_t> trace;
  std::vector<uint16_t> dtrace1;
  std::vector<uint16_t> dtrace2;
  if(TreeEna){
    fTree->Branch("trace"   , &trace       );
    fTree->Branch("dprobe1" , &dtrace1     );
    fTree->Branch("dprobe2" , &dtrace2     );
    fTree->Branch("tstamp" , &tstamp_10ns, "tstamp/l" );
    fTree->Branch("cfd"    , &cfd        , "cfd/F"    );
    fTree->Branch("board"  , &boardId    , "board/s"  );
    fTree->Branch("channel", &channelId  , "channel/s");
    fTree->Branch("qshort" , &qshort     , "qshort/s" );
    fTree->Branch("qslong" , &qlong      , "qlong/s" );
    fTree->Branch("time"   , &time       , "time/s"   );
  }
  
  // RAW DATA file has a header containing the number of boards and the
  // AGAVA timestamp corresponding to the T=0 of the CAEN boards 
  // i.e. the start of the run
  uint64_t tsOffset;
  uint16_t nbOfBoards;
  uint32_t *sizeFrame = new uint32_t;
  fread(sizeFrame,1,sizeof(uint32_t),fInput);
  *sizeFrame = *sizeFrame & 0x0FFFFFFF;
  uint32_t countingRate[6][16];
  for(uint32_t bd = 0 ; bd < 6 ; bd++){
    for(uint32_t ch = 0 ; ch < 16 ; ch++){
      countingRate[bd][ch]=0;
    }
  }
  //  std::cout << "Size to read = " << *sizeFrame << std::endl;
  char *buffer = new char[512*512];

  // READING THE XDAQ HEADER 
  fread(buffer,*sizeFrame-1,sizeof(uint32_t),fInput);
  
  uint32_t *convBuff= (uint32_t*) buffer;
  uint64_t fGlobalOffset = (uint64_t)(convBuff[2])<<32|convBuff[3];
  
  printf("fGlobalOffset = %llx\n",fGlobalOffset);
  time_t rawtime = convBuff[4];
  struct tm  ts;
  char       buf[80];

  // Format time, "ddd yyyy-mm-dd hh:mm:ss zzz"
  ts = *localtime(&rawtime);
  strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", &ts);
  printf("%s\n", buf);
//  for(uint32_t j = 0 ; j < *sizeFrame-1 ; ++j)
//    printf("0x%08X\n",convBuff[j]);
  uint32_t aggregateSize =0;
  uint32_t nbAggr=0;
  uint32_t offset = 5;
  uint32_t channelAggrSize = 0;
  FILE *pFile;
  pFile = fopen(argv[2],"wb");
  
  while((nbMaxAggr==0 && !feof(fInput)) || (nbMaxAggr > 0 && nbAggr < nbMaxAggr) ){
    fread(sizeFrame,1,sizeof(uint32_t),fInput);
    aggregateSize = (*sizeFrame&0x0FFFFFFF)-1;
    //    printf("Size to read = %u (0x%08X)\n",aggregateSize,*sizeFrame);
    fread(buffer,aggregateSize,sizeof(uint32_t),fInput);

    //decoding the aggregate: 
    convBuff = (uint32_t*) buffer;
    boardId = (convBuff[0]&0xF8000000)>>27;
    std::bitset<8> channelMask = convBuff[0]&0xFF;
    // we need to reconstruct the channel number from the couple id for the 
    // 16 channels boards  
    uint8_t couples[8]={0,0,0,0,0,0,0,0};
    int index = 0 ;
    std::bitset<16> extraBit;
    for(uint16_t bit = 0 ; bit < 8 ; bit++){ 
      if(channelMask.test(bit)){
	//	std::cout << "Couple " << bit << " activated" << std::endl;
	couples[index]=bit;
	index++;
      }
    }
    uint64_t tstamp_10ns  = 0;
    uint32_t fineTS       = 0;
    uint32_t extras,extras2;
    uint32_t coupleAggregateSize = 0 ;
    uint32_t startingPos  = 3;
    uint32_t fNsPerTStamp = 4;
	int	tsLocalOffset =0;
//        std::cout << "Board Id = " << boardId << std::endl;
    // passing to channel aggregates
    for(uint8_t cpl = 0 ; cpl < channelMask.count() ; cpl++)
    {
        coupleAggregateSize = convBuff[startingPos]&0x3FFFFF;
        dataForm.SetDataFormat(convBuff[startingPos+1]);
        //dataForm.printDataFormat();
        for(uint16_t evt = 0 ; evt < (coupleAggregateSize-2)/dataForm.evtSize();evt++)
        {
            cfd = 0.;
            tstamp  = static_cast<uint64_t>((convBuff[startingPos+2+dataForm.evtSize()*evt])&0x7FFFFFFF);
            if(dataForm.enaExtras())
            {
                extras  = static_cast<uint16_t>(((convBuff[startingPos+dataForm.evtSize()-1+dataForm.evtSize()*evt+2])&0x1F0000)>>16);
                extraBit = extras;
                isat = extraBit[ 4];
                tsat = extraBit[10];
                pur  = extraBit[ 9];
                extras2 = convBuff[startingPos+dataForm.evtSize()-2+dataForm.evtSize()*evt+2];
                switch(dataForm.confExtras()) 
                {
                    case 0: // Extended Time Stamp [31:16] ; baseline*4 [15:0]
                        tstamp = ((uint64_t)(extras2&0xFFFF0000))<<15 | (uint64_t)tstamp;
                        break;
                    case 1: // Extended Time stamp [31:16] ; flags [15:0]
                        tstamp = ((uint64_t)(extras2&0xFFFF0000))<<15 | (uint64_t)tstamp;
                        break;
                    case 2: // Extended Time stamp [31:16] ; Flags [15:10] ; Fine Time Stamp [9:0]
                        fineTS = static_cast<uint16_t>(extras2&0x3FF);
                        cfd = fNsPerTStamp/1.024* fineTS; // to have it in ps units
                //	        std::cout << fineTS << " " << cfd << std::endl;
                        tstamp = ((uint64_t)(extras2&0xFFFF0000))<<15 | (uint64_t)tstamp;
                        break;
                    case 4: // Lost Trigger Counter [31:16] ; Total Trigger [15:0]
                        break;
                    case 5: // Positive zero crossing [31:16] ; Negative zero crossing [15:0]
                        break;
                    case 7: // Debug fixed value;
                        break;
                    default:
                        break;
                }
            }
            //      statsData.time[cpl] = tstamp*2/pow( 10, 9 );   // must be in ns
            channelId          = static_cast<uint32_t>(((convBuff[startingPos+2+dataForm.evtSize()*evt]&0x80000000)>>31)
                                +couples[cpl]*2);
            qshort = static_cast<uint16_t>((convBuff[startingPos+dataForm.evtSize()-1+dataForm.evtSize()*evt+2])&0x7FFF);
            qlong  = static_cast<uint16_t>((convBuff[startingPos+dataForm.evtSize()-1+dataForm.evtSize()*evt+2])>>16);
            uint16_t framePos = startingPos+3+dataForm.evtSize()*evt;
            trace.clear();
            trace.resize(dataForm.numSamples());
            dtrace1.clear();
            dtrace1.resize(dataForm.numSamples());
            dtrace2.clear();
            dtrace2.resize(dataForm.numSamples());
            for (uint16_t samp=0 ; samp<dataForm.numSamples()/2 ; samp++){
            trace[2*samp]     = static_cast<uint16_t>(convBuff[framePos]&0x3FFF);
            trace[2*samp+1]   = static_cast<uint16_t>((convBuff[framePos]>>16)&0x3FFF);
            dtrace1[2*samp]   = static_cast<bool>((convBuff[framePos]&0x8000)>>15);
            dtrace1[2*samp+1] = static_cast<bool>((convBuff[framePos]&0x80000000)>>31);
            dtrace2[2*samp]   = static_cast<bool>((convBuff[framePos]&0x4000)>>14);
            dtrace2[2*samp+1] = static_cast<bool>((convBuff[framePos]&0x40000000)>>30);
            //printf("%3i %08X --> %04X %04X\n",samp,convBuff[framePos],trace[samp],trace[samp+1]);
            framePos++;
            }
            if(debug){
            std::cout << "--------------------------\n"; 
            int line=0;
            std::cout << "--------------------------\nL"<<line<<'\t'; 
            for (uint16_t samp=0 ; samp<dataForm.numSamples() ; samp++){
                if(dtrace2[samp] && !dtrace1[samp])
                printf(ANSI_COLOR_GREEN "%5i   "ANSI_COLOR_RESET,trace[samp]);
                else if(!dtrace2[samp] && dtrace1[samp])
                printf(ANSI_COLOR_RED "%5i   "ANSI_COLOR_RESET,trace[samp]);
                else if (dtrace2[samp] && dtrace1[samp])
                printf(ANSI_COLOR_MAGENTA "%5i   "ANSI_COLOR_RESET,trace[samp]);
                else
                printf("%5i   ",trace[samp]);
                if((samp+1)%16==0) std::cout << std::endl << "L"<<++line<<'\t';
            }
            std::cout << "--------------------------\n"; 
            }
            //printf("%3d %6d %6d\n",channelId,qshort,qlong);
            //	dbug_line;
            if(boardId == 1 && (channelId%4==0 || (channelId-1)%4==0)) { /*std::cout << boardId << " " << channelId << std::endl;*/ tsLocalOffset = -60;}
            else tsLocalOffset =0;
            tstamp_10ns        = tstamp * fNsPerTStamp / 10 ;
            cfd               += (tstamp*fNsPerTStamp - tstamp_10ns*10)*1000;
            tstamp_10ns       += fGlobalOffset + tsLocalOffset;

            energy  = static_cast<uint16_t>((convBuff[startingPos+dataForm.evtSize()-1+dataForm.evtSize()*evt+2])&0x7FFF);
        //	if(energy>700)std::cout << " extras = " << std::hex << extras 
        //				<< " energy = " << std::dec << energy << std::endl;
            
            time    = static_cast<uint16_t>(std::floor(cfd));

            if(energy>0)countingRate[boardId][channelId]++;
            fTree->Fill();
            trace.clear();
            dtrace1.clear();
            dtrace2.clear();
        }
        
        startingPos+=coupleAggregateSize;
    }
    
    //     for(uint32_t j = 0 ; j < aggregateSize+10 ; j++)
    //       printf("\t %08X\n",convBuff[j]);
    memset(buffer,512*512,0); //this is ineffective
    nbAggr++;
        //std::cout << nbAggr << std::endl;
  }
  
  std::cout << "Number of aggregate = " << nbAggr << std::endl;
  std::cout << "Number of events with Energy > 0 channels: "<< std::endl;
  for(uint32_t bd = 0 ; bd < 6 ; bd++){
    for(uint32_t ch = 0 ; ch < 16 ; ch++){
      if(countingRate[bd][ch]>0){
	printf("Board %2i Channel %2i NbEvents = %8i\n",bd,ch,countingRate[bd][ch]);
      }
    }
  }
  //  fTree->Write(0,TObject::kOverwrite);
  fOutput->Write();
  fOutput->Close();
  fclose(fInput);
  fclose(pFile);
//  delete [] buffer;
//  delete [] convBuff;
//  delete sizeFrame;
  return 0;
  }
