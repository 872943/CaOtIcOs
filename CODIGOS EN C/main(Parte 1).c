#include "funciones.h"


int main() {
    double dt=0.1;
    int N_pasos=50;
    int N_redes=100;
    int rede_ini=620;




    double betta_inicial=0;
    double betta_final=3;
    double K_inicial=0;
    double K_final=10;
    double delta_K=0.1;
    double delta_betta=0.05;
    
    double betta_actual=betta_inicial;
    double K_actual=K_inicial;
    
    while(K_actual<K_final){
        printf("+1");
        while(betta_actual<betta_final){
            frac_polarizado(N_redes,rede_ini, K_actual, betta_actual, "RESULTADOS_MAPA", N_pasos,dt);
            betta_actual=betta_actual+delta_betta;
        }
        K_actual=K_actual+delta_K;
        printf("+1");
    }
    

}