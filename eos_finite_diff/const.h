//allgemein
double const c_g = 6.6743 * pow(10, -11);  // Gravitationskonstante G
double const c_c = 299792458.0;              // Lichtgeschwindigkeit c
double const c_gamma = (1+1.0/0.5);          // Der Exponent des Polytropen EOS mit dem Polytropenindex n, rho = (P/K)^gamma    
double const c_k = 2.*pow(10,19);             // Konstante K 2*pow(10,8) für n=1;
double const c_pi= 3.14159265359;            // pi
double hbar = 1.054571817*pow(10,-34);        //h quer
double echarge = 1.602176634*pow(10,-19);     //elementarladung
double Nc = 3.;                               // Anzahl der Farbladungen 
double nu_si=(pow(10,36)*pow(hbar*c_c,-3)*pow(echarge,4)); //Faktor für Umrechunung NU zu SI beim Druck
double si_mevfm3 = pow(10, -51) / echarge;
double si_geom = pow(c_c, -4) * c_g;