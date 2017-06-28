#include <stdlib.h>
#include <iostream>
#include <time.h>
double RandomNumber(int L)
{
	srand(0);
	double a = rand()%200;
	return a/200;
}
int main()
{
	srand(time(NULL));
	double R[10];
	for (int i =0; i<5; i++)
	{
		R[i]=rand() % 100;
	std::cout << R[i]/100 << "\n";
	}
std:: cout <<"\n" << RandomNumber(2);
}
