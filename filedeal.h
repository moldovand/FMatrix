#ifndef _filedeal_h_
#define _filedeal_h_

//�ե�����ݥ��󥿤򼡤ιԤ���Ƭ���֤�
void NextLine(FILE* fp);

//�ե�����ݥ��󥿤򼡤ο�����Ƭ���֤�
void NextNumber(FILE* fp);

//�ݥ���ȥǡ�����٥��ȥ�˳�Ǽ
void Text2(char* fname, double**& vec, const int& pnum, 
	   double& xm, double& ym, double& Lx, double& Ly);

#endif //_filedeal_h_
