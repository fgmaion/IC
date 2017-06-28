double Kernel(double q, double h)
{
	double R = 0;
	if (0<=q && q<1)
		{
			R = (pow((2-q),3)-4*pow((1-q),3))/(6*h);
		}
	if ((1<=q && q<2 ))
		{
			R = (pow((2-q),3))/(6*h);
		}
	if (2<=q)
		{
			R = 0;
		}
	return (R);
}
