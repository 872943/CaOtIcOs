#include "funciones.h"


int main() {
    double dt=0.1;
    int N_pasos=50;
    int N_redes=100;
    int rede_ini=620;




    double betta_inicial=0;
    double betta_final=2;
    double K_inicial=0;
    double K_final=20;
    double delta_K=0.2;
    double delta_betta=0.05;
    
    double betta_actual=betta_inicial;
    double K_actual=K_final;
    
    while(K_actual>K_inicial){
        printf("+1");
        betta_actual=betta_inicial;
        while(betta_actual<betta_final){
            frac_polarizado(N_redes,rede_ini, K_actual, betta_actual, "RESULTADOS_MAPA_Irene", N_pasos,dt);
            betta_actual=betta_actual+delta_betta;
        }
        K_actual=K_actual-delta_K;
        printf("+1");
    }
    

}