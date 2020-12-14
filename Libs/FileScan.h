#ifndef _FileScan_h_
#define _FileScan_h_

#define IDETA  1
#define DIGIT  2
#define P_NUM  3

int FeedComment(FILE *fp);
int ScanForScan(FILE *fp);

#endif  // _FileScan_h_
