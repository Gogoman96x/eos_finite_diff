using namespace std;
int b=-1;
//parameter  für k
double Lambda = 1.; //ev obere Grenze k
//double kir = 0.1;  //ev unter Grenze k
//double delk = 0.0001; //ev Schrittweite von k
//double k=Lambda;
int a=0;
double k_check;



//mfa
//double a1=pow(901.09e-3,2);
//double a2=-2*5.38;

//FRG parameter
double a1=pow(711.98e-3,2);
double a2 = 20.22;
double hl=pow(120.73e-3,3);
double yk_g=6.5;
double mu,check,T;



typedef vector< double> state_type;
typedef vector<vector< double>> state_type2;


//Globale Variablen für Interpol
double xmin=pow(1.e-3,1);
double xmax=pow(170.e-3,1);
double xn=xmin;
int const points=40;
double xarr[points];
double xarr2[points];
double delxarr[points];
double delxarr2[points];
int const n_th=points-2;
double c_th[points-1];
double b_th[points];
double a_th[points-1];
double derivs[2][points];
double derivs2[2][points];
double derivs3[2][points];
double finite_diff_num[5];
double finite_diff_dem[5];
double finite_diff[5];
double finite_diff2_num[5];
double finite_diff2_dem[5];
double finite_diff2[5];
double algo1[points - 1]; //c prime vom thomas algorithmus
double algo2[points]; // d prime vom Thomas algorithmus
double algo3[points]; //mi berechnete
double abl1[points]; // ersten Ableitungen am Punkt xi
double abl2[points]; // zweite Ableitungen am Punkt xi
double funcdata[points]; // Werte dj hier rein
double data1[points];
double data2[points];
state_type datatest;



void funcinterp(state_type datay,double datax[], int j);
double func2(double x1, double x2, double x3,double x4,double x5);
double func1(double x1);
double testfunc(double x);
/*generiert die Punkte yi, da über k integriert wird,
  bleibt die Funktion nicht konstant. x2 ist Term der durch die flow eq. folgt.*/
double func2(double x1, double x2, double x3,double x4,double x5){
    double F2=((x3-x2)/x5)-((x2-x1)/x4);
    return F2;
}
//alle arrays werden hier gefüllt
void def_var(){
    xarr[0]=xmin;
    
    double delx=(xmax-xmin)/(points-1);
    for(int i=1;i<points;i++){
        xarr[i]=xmin+i*delx;
        xarr2[i] = pow(xarr[i], 2);
    }

    for(int i=0;i<points;i++){
        datatest.push_back(testfunc(xarr2[i]));
        //data1[i]=func1(xn,0);
        data2[i]=func1(xarr2[i]);//hier
        //floweq[i]=0;
        //xn+=delxarr[i];
    }

    //funcinterp(data2,xarr);
}


//void funktion berechnet die zweiten und ersten Ableitungen. 
void funcinterp(state_type datay_o,double datax[],int j){
    state_type datay;

    if (j == 1)
    {
        for (int i = 0; i < points; i++)
        {
            datay.push_back(datay_o[i]);
        }
    }

    if (j == 2)
    {
        for (int i = points; i < 2*points; i++)
        {
            datay.push_back(datay_o[i]);
        }
    }

    if (j == 3)
    {
        for (int i = 2*points; i < 3*points; i++)
        {
            datay.push_back(datay_o[i]);
        }
    }

    for(int i=0;i<points-1;i++){
        delxarr2[i]=datax[i+1]-datax[i];
    }
    b_th[0]=delxarr2[0]/3.;
    b_th[points-1]=delxarr2[points-2]/3.;
    for(int i=1;i<points-1;i++){
        b_th[i]=(delxarr2[i-1]+delxarr2[i])/3.;
    }
    c_th[0]=delxarr2[0]/6.;
    a_th[points-2]=delxarr2[points-2]/6.;
    for(int i=0;i<points-2;i++){
        a_th[i]=delxarr2[i]/6;
        c_th[i+1]=delxarr2[i+1]/6;
    }

    

    //abl1[0] = (datay[1] - datay[0]) / delxarr2[0];
    //abl1[points-1]= (datay[points-1] - datay[points-2]) / delxarr2[points-2];


    abl1[0] = (2 * datax[0] - datax[1] - datax[2]) / ((datax[1] - datax[0]) * (datax[2] - datax[0])) * datay[0]+
    (datax[2] - datax[0]) / ((datax[1] - datax[0]) * (datax[2] - datax[1])) * datay[1]+
    (datax[0] - datax[1]) / ((datax[2] - datax[1]) * (datax[2] - datax[0])) * datay[2];
    
    abl1[points - 1] = (datax[points - 1] - datax[points - 2]) / ((datax[points - 1] - datax[points - 3]) * (datax[points - 2] - datax[points - 3])) * datay[points - 3] +
        (datax[points - 3] - datax[points - 1]) / ((datax[points - 2] - datax[points - 3]) * (datax[points - 1] - datax[points - 2])) * datay[points - 2] +
        (2 * datax[points - 1] - datax[points - 2] - datax[points - 3]) / ((datax[points - 1] - datax[points - 2]) * (datax[points - 1] - datax[points - 3])) * datay[points - 1];



    for(int i=0;i<points;i++){
        algo2[i]=0;
        algo3[i]=0;
    }
    for(int i=0;i<points-1;i++){
        algo1[i]=0;
    }




    funcdata[0]=(datay[1]-datay[0])/delxarr2[0]-abl1[0];
    funcdata[points-1]=-(datay[points-1]-datay[points-2])/delxarr2[points-2]+abl1[points-1];

    for(int i=1;i<points-1;i++){//Werte d_j werden in Array geschrieben.
        funcdata[i]=func2(datay[i-1],datay[i], datay[i+1],delxarr2[i-1],delxarr2[i]);
    }

    algo1[0]=c_th[0]/b_th[0];
    for(int i=1;i<points-1;i++){
        algo1[i]=c_th[i]/(b_th[i]-algo1[i-1]*a_th[i-1]);
    }

   
    algo2[0]=funcdata[0]/b_th[0];
    for(int i=1;i<points;i++){
            algo2[i]=(funcdata[i]-a_th[i-1]*algo2[i-1])/(b_th[i]-a_th[i-1]*algo1[i-1]);
    }

    algo3[points-1]=algo2[points-1];

    for(int i=points-2; i>=0;i--){
        algo3[i]=algo2[i]-algo1[i]*algo3[i+1];
    }


    for(int i=0;i<points;i++){
        abl2[i]=algo3[i];
    }

    for(int i=1;i<points-1;i++){
        abl1[i]=(datay[i+1]-datay[i])/delxarr2[i] -delxarr2[i]*(abl2[i+1]+2*abl2[i])/6;
    }
    if (j == 1)
    {
        for (int i = 0; i < points; i++) {
            derivs[0][i] = abl1[i];
            derivs[1][i] = abl2[i];
        }
    }

    if (j == 2)
    {
        for (int i = 0; i < points; i++) {
            derivs2[0][i] = abl1[i];
            derivs2[1][i] = abl2[i];
        }
    }

    if (j == 3)
    {
        for (int i = 0; i < points; i++) {
            derivs3[0][i] = abl1[i];
            derivs3[1][i] = abl2[i];
        }
    }

}



