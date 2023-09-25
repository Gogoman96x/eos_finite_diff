#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>
#include <fstream>
#include "interpol_newfrg2.h"
#include "const.h"
#include "integr.h"
#include "spline.h"


using namespace std;
using namespace boost::numeric::odeint;



double func1(double x1){
        double F1=1./2*a1*x1+a2/8.*pow(x1,2);
        return F1;
}


double testfunc(double x){
    double test=pow(x,2);
    return test;
}

state_type xtest;

double funcflow(double k, double mu, double sig, double t, state_type x, int i)
{   

    double ef = sqrt(pow(k, 2) + pow(yk_g / 2 * sig, 2));
    double eb1 = sqrt(pow(k, 2) + 2 * derivs[0][i] + 4 * pow(sig, 2) * derivs[1][i]);
    double eb2 = sqrt(pow(k, 2) + 2 * derivs[0][i]);

    double F1 = pow(k, 4) / (12 * pow(c_pi, 2)) *
        (1 / eb1 * 1 / tanh(eb1 / (2 * t)) + 3 / eb2 * 1 / tanh(eb2 / (2 * t)));
    
    double F2 = -pow(k, 4) / (3 * pow(c_pi, 2)) *
        (Nc * (1 / ef) * (tanh((ef - mu) / (2 * t)) + tanh((ef + mu) / (2 * t))));
    double F =F1+F2;
    return F;
}

double funcflow_vac(double k, double mu, double sig, double t, state_type x, int i)
{

    double ef = sqrt(pow(k, 2) + pow(yk_g / 2 * sig, 2));
    double eb1 = sqrt(pow(k, 2) + 2 * derivs[0][i] + 4 * pow(sig, 2) * derivs[1][i]);
    double eb2 = sqrt(pow(k, 2) + 2 * derivs[0][i]);
    double F1_t0;
    double F2_t0;

    F1_t0 = pow(k, 4) / (12 * pow(c_pi, 2)) * (1 / eb1 + 3 / eb2);
    if (ef >= mu) {
        F2_t0 = pow(k, 4) / (12 * pow(c_pi, 2)) * (-8 * Nc * (1 / ef));
    }
    else if (ef < mu)
    {
        F2_t0 = 0;
    }
    double F = F1_t0 + F2_t0;
    return F;
}


void floweq(const state_type& x, state_type& dxdt, const double k)
{
    funcinterp(x, xarr2,1);

    for (int i = 0; i <points; i++)
    {
        dxdt[i] = funcflow(k, mu, xarr[i], T, x, i);
    }
}

void floweq_vac(const state_type& x, state_type& dxdt, const double k)
{
    funcinterp(x, xarr2, 1);
    for (int i = 0; i < points; i++)
    {
        dxdt[i] = funcflow_vac(k, mu, xarr[i], T, x, i);
    }
}

void write_cout(const state_type& x, const double t)
{
    cout << t << '\t';
        cout << x[0] << '\t';
    
    cout << endl;
}



typedef runge_kutta_dopri5< state_type > stepper_type;
typedef euler<state_type> eu;
int main()
{
    double mu_og = 2e-4;
    double T_og = 0.001;
    double omeg_zero;
    ofstream file1;
    file1.open("frg_t1_k80.dat");
    def_var();
    //  U(0,0) wird hier berechnet
    state_type x_vac;
    for (int i = 0; i < points; i++)
    {
        x_vac.push_back(data2[i]);
    }
    mu = 0;
    integrate_adaptive(make_controlled(1e-14, 1e-14, stepper_type()),
        floweq_vac, x_vac, 1., 0.08, -0.00001);


    while (mu <= 0.5)
    {
        //Es werden vier Punkte berechnet U(T0+delT,mu0),U(T0-delT,mu0)...
        T = T_og;
        mu = mu_og;
        state_type x;
        state_type x_pdelt;
        state_type x_mdelt;
        state_type x_pdelmu;
        state_type x_mdelmu;

        //fünf Anfangsbedingungen werden hier deklariert
        for (int i = 0; i < points; i++)
        {
            x.push_back(data2[i]);
            x_pdelt.push_back(data2[i]);
            x_mdelt.push_back(data2[i]);
            x_pdelmu.push_back(data2[i]);
            x_mdelmu.push_back(data2[i]);

        }
        //Die fünf PDE's werden hier gelöst. Die Parameter Werte T und mu werden jeweils variiert.
        //############################################################################################
        integrate_adaptive(make_controlled(1e-14, 1e-14, stepper_type()),
            floweq, x, 1., 0.08, -0.00001);
        T = T_og+2e-4;
        integrate_adaptive(make_controlled(1e-14, 1e-14, stepper_type()),
            floweq, x_pdelt, 1., 0.08, -0.00001);
        T = T_og-2e-4;
        
        integrate_adaptive(make_controlled(1e-14, 1e-14, stepper_type()),
            floweq, x_mdelt, 1., 0.08, -0.00001);
        T = T_og;
        mu = mu_og + 2e-4;
        integrate_adaptive(make_controlled(1e-14, 1e-14, stepper_type()),
            floweq, x_pdelmu, 1., 0.08, -0.00001);
        mu = mu_og - 2e-4;
        integrate_adaptive(make_controlled(1e-14, 1e-14, stepper_type()),
            floweq, x_mdelmu, 1., 0.08, -0.00001);
        //############################################################################################

        //Der symmetrie breaking term muss mit betrachtet werden zur Bestimmung von Omega(T,mu)
        state_type x_fin;
        state_type y_fin;
        state_type y_pdelt;
        state_type y_mdelt;
        state_type y_pdelmu;
        state_type y_mdelmu;
        state_type y_vac;

        //Vektoren der Koordinaten für Splines werden gefüllt.
        for (int i = 0; i < points; i++)
        {
            y_vac.push_back(x_vac[i] - hl * xarr[i]);
            y_fin.push_back(x[i] - hl * xarr[i]);
            y_pdelt.push_back(x_pdelt[i] - hl * xarr[i]);
            y_mdelt.push_back(x_mdelt[i] - hl * xarr[i]);
            y_pdelmu.push_back(x_pdelmu[i] - hl * xarr[i]);
            y_mdelmu.push_back(x_mdelmu[i] - hl * xarr[i]);
            x_fin.push_back(xarr[i]);
        }
        //splines werden definiert
        tk::spline s0(x_fin, y_vac);
        tk::spline s(x_fin, y_fin);
        tk::spline s1(x_fin, y_pdelt);
        tk::spline s2(x_fin, y_mdelt);
        tk::spline s3(x_fin, y_pdelmu);
        tk::spline s4(x_fin, y_mdelmu);


        //Parameter zur Bestimmung der Minima vom Spline werden definiert
        double x_spline = xmin;//stepper
        double xmin0_spline = xmin;
        double xmin_spline = xmin;
        double xmin1_spline = xmin;
        double xmin2_spline = xmin;
        double xmin3_spline = xmin;
        double xmin4_spline = xmin;
        double ymin0_spline = s0(x_spline);
        double ymin_spline = s(x_spline);
        double ymin1_spline = s1(x_spline);
        double ymin2_spline = s2(x_spline);
        double ymin3_spline = s3(x_spline);
        double ymin4_spline = s4(x_spline);

        //Minimum wird berechnet von den Splines 
        while (x_spline < xmax)
        {
            if (s0(x_spline) < ymin0_spline)
            {
                ymin0_spline = s0(x_spline);
                xmin0_spline = x_spline;
            }
            if (s(x_spline) < ymin_spline)
            {
                ymin_spline = s(x_spline);
                xmin_spline = x_spline;
            }
            if (s1(x_spline) < ymin1_spline)
            {
                ymin1_spline = s1(x_spline);
                xmin1_spline = x_spline;
            }
            if (s2(x_spline) < ymin2_spline)
            {
                ymin2_spline = s2(x_spline);
                xmin2_spline = x_spline;
            }
            if (s3(x_spline) < ymin3_spline)
            {
                ymin3_spline = s3(x_spline);
                xmin3_spline = x_spline;
            }
            if (s4(x_spline) < ymin4_spline)
            {
                ymin4_spline = s4(x_spline);
                xmin4_spline = x_spline;
            }
            x_spline += 0.00001;
        }

        omeg_zero = ymin0_spline;
        double p = omeg_zero - s(xmin_spline);
        double delt = (ymin1_spline - ymin2_spline) / (2 * 2e-4);
        double delmu = (ymin3_spline - ymin4_spline) / (2 * 2e-4);
        
        file1 <<mu_og*1000<<'\t'<<xmin_spline*1000<<
            '\t' << p * nu_si * si_mevfm3 << '\t' << (-p - T_og * delt - mu_og * delmu) * nu_si * si_mevfm3<<endl;
        
        cout << mu_og << endl;
        mu_og += 0.001;
    }
    file1.close();
    
    return 0;
}