void metric_JP(double spin, double epsilon, 
			   double r, double theta, double mn[4][4]);

void Christoffel(double spin, double epsilon, double w1, 
	             double w2, double CS[][4][4])
{
	double spin2 = spin*spin;
	double g[4][4];
	double gg;
	double gr[4][4];
//	double grr[4][4];
	double gl[4][4];
//	double gll[4][4];
	double dw1, dw2, temp;
	double invg[4][4];
	double Dg[4][4][4];
	
	/* ----- set dw1 and dw2 ----- */
	dw1  = 0.01*w1;
	temp = w1 + dw1;
	dw1  = temp - w1;
	dw2  = 0.01;
	temp = w2 + dw2;
	dw2  = temp - w2;

	/* ----- compute inverse metric ----- */
	metric_JP(spin, epsilon, w1, w2, g);
	gg = g[0][0]*g[3][3] - g[0][3]*g[0][3];
	
	invg[0][0] = g[3][3]/gg; 
	invg[0][3] = - g[0][3]/gg;
	invg[1][1] = 1/g[1][1];
	invg[2][2] = 1/g[2][2];
	invg[3][0] = invg[0][3]; 
	invg[3][3] = g[0][0]/gg;
	
	/* ----- compute metric derivatives ----- */
	metric_JP(spin, epsilon, w1 + dw1, w2, gr);
	metric_JP(spin, epsilon, w1 - dw1, w2, gl);
	
	Dg[0][0][1] = 0.5*(gr[0][0] - gl[0][0])/dw1;
	Dg[0][3][1] = 0.5*(gr[0][3] - gl[0][3])/dw1;
	Dg[1][1][1] = 0.5*(gr[1][1] - gl[1][1])/dw1;
	Dg[2][2][1] = 0.5*(gr[2][2] - gl[2][2])/dw1;
	Dg[3][0][1] = Dg[0][3][1];
	Dg[3][3][1] = 0.5*(gr[3][3] - gl[3][3])/dw1;
	
	metric_JP(spin, epsilon, w1, w2 + dw2, gr);
	metric_JP(spin, epsilon, w1, w2 - dw2, gl);
	
	Dg[0][0][2] = 0.5*(gr[0][0] - gl[0][0])/dw2;
	Dg[0][3][2] = 0.5*(gr[0][3] - gl[0][3])/dw2;
	Dg[1][1][2] = 0.5*(gr[1][1] - gl[1][1])/dw2;
	Dg[2][2][2] = 0.5*(gr[2][2] - gl[2][2])/dw2;
	Dg[3][0][2] = Dg[0][3][2];
	Dg[3][3][2] = 0.5*(gr[3][3] - gl[3][3])/dw2;
	/*
	metric(w1+dw1,w2,gr);
	metric(w1+dw1+dw1,w2,grr);
	metric(w1-dw1,w2,gl);
	metric(w1-dw1-dw1,w2,gll);
	
	Dg[0][0][1] = (8*gr[0][0] - grr[0][0] - 8*gl[0][0] + gll[0][0])/dw1/12;
	Dg[0][3][1] = (8*gr[0][3] - grr[0][3] - 8*gl[0][3] + gll[0][3])/dw1/12;
	Dg[1][1][1] = (8*gr[1][1] - grr[1][1] - 8*gl[1][1] + gll[1][1])/dw1/12;
	Dg[2][2][1] = (8*gr[2][2] - grr[2][2] - 8*gl[2][2] + gll[2][2])/dw1/12;
	Dg[3][0][1] = Dg[0][3][1];
	Dg[3][3][1] = (8*gr[3][3] - grr[3][3] - 8*gl[3][3] + gll[3][3])/dw1/12;
	
	metric(w1,w2+dw2,gr);
	metric(w1,w2+dw2+dw2,grr);
	metric(w1,w2-dw2,gl);
	metric(w1,w2-dw2-dw2,gll);
	
	Dg[0][0][2] = (8*gr[0][0] - grr[0][0] - 8*gl[0][0] + gll[0][0])/dw2/12;
	Dg[0][3][2] = (8*gr[0][3] - grr[0][3] - 8*gl[0][3] + gll[0][3])/dw2/12;
	Dg[1][1][2] = (8*gr[1][1] - grr[1][1] - 8*gl[1][1] + gll[1][1])/dw2/12;
	Dg[2][2][2] = (8*gr[2][2] - grr[2][2] - 8*gl[2][2] + gll[2][2])/dw2/12;
	Dg[3][0][2] = Dg[0][3][2];
	Dg[3][3][2] = (8*gr[3][3] - grr[3][3] - 8*gl[3][3] + gll[3][3])/dw2/12;
	*/
	
	/* ----- compute Christoffel symbols ----- */
	CS[0][0][0] = 0;
	CS[0][0][1] = invg[0][0]*Dg[0][0][1] + invg[0][3]*Dg[0][3][1];
	CS[0][0][2] = invg[0][0]*Dg[0][0][2] + invg[0][3]*Dg[0][3][2];
	CS[0][0][3] = 0;
	CS[0][1][0] = CS[0][0][1];
	CS[0][1][1] = 0;
	CS[0][1][2] = 0;
	CS[0][1][3] = invg[0][0]*Dg[0][3][1] + invg[0][3]*Dg[3][3][1];
	CS[0][2][0] = CS[0][0][2];
	CS[0][2][1] = 0;
	CS[0][2][2] = 0;
	CS[0][2][3] = invg[0][0]*Dg[0][3][2] + invg[0][3]*Dg[3][3][2];
	CS[0][3][0] = 0;
	CS[0][3][1] = CS[0][1][3];
	CS[0][3][2] = CS[0][2][3];
	CS[0][3][3] = 0;
	
	CS[1][0][0] = - invg[1][1]*Dg[0][0][1];
	CS[1][0][1] = 0;
	CS[1][0][2] = 0;
	CS[1][0][3] = - invg[1][1]*Dg[0][3][1];
	CS[1][1][0] = 0;
	CS[1][1][1] = invg[1][1]*Dg[1][1][1];
	CS[1][1][2] = invg[1][1]*Dg[1][1][2];
	CS[1][1][3] = 0;
	CS[1][2][0] = 0;
	CS[1][2][1] = CS[1][1][2];
	CS[1][2][2] = - invg[1][1]*Dg[2][2][1];
	CS[1][2][3] = 0;
	CS[1][3][0] = CS[1][0][3];
	CS[1][3][1] = 0;
	CS[1][3][2] = 0;
	CS[1][3][3] = - invg[1][1]*Dg[3][3][1];
	
	CS[2][0][0] = - invg[2][2]*Dg[0][0][2];
	CS[2][0][1] = 0;
	CS[2][0][2] = 0;
	CS[2][0][3] = - invg[2][2]*Dg[0][3][2];
	CS[2][1][0] = 0;
	CS[2][1][1] = - invg[2][2]*Dg[1][1][2];
	CS[2][1][2] = invg[2][2]*Dg[2][2][1];
	CS[2][1][3] = 0;
	CS[2][2][0] = 0;
	CS[2][2][1] = CS[2][1][2];
	CS[2][2][2] = invg[2][2]*Dg[2][2][2];
	CS[2][2][3] = 0;
	CS[2][3][0] = CS[2][0][3];
	CS[2][3][1] = 0;
	CS[2][3][2] = 0;
	CS[2][3][3] = - invg[2][2]*Dg[3][3][2];
	
	CS[3][0][0] = 0;
	CS[3][0][1] = invg[3][3]*Dg[0][3][1] + invg[3][0]*Dg[0][0][1];
	CS[3][0][2] = invg[3][3]*Dg[0][3][2] + invg[3][0]*Dg[0][0][2];
	CS[3][0][3] = 0;
	CS[3][1][0] = CS[3][0][1];
	CS[3][1][1] = 0;
	CS[3][1][2] = 0;
	CS[3][1][3] = invg[3][3]*Dg[3][3][1] + invg[3][0]*Dg[0][3][1];
	CS[3][2][0] = CS[3][0][2];
	CS[3][2][1] = 0;
	CS[3][2][2] = 0;
	CS[3][2][3] = invg[3][3]*Dg[3][3][2] + invg[3][0]*Dg[0][3][2];
	CS[3][3][0] = 0;
	CS[3][3][1] = CS[3][1][3];
	CS[3][3][2] = CS[3][2][3];
	CS[3][3][3] = 0;

	return ;
}
