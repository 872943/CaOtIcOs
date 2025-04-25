#include "funciones.h"


int main() {

    srand(time(NULL)); // Inicializa aleatoriedad
    double dt=0.1;
    int N_pasos=500;
    double K=10;;
    int N_redes=10;
    int rede_ini=1;

    double betta_actual=-0.1;
    double betta_final=2;
    double delta_betta=0.1;

    while(betta_actual<betta_final){
        betta_actual=betta_actual+delta_betta;
        frac_polarizado(N_redes,rede_ini,K,betta_actual,N_pasos,dt);
        calcular_fraccion_polarizados(K,betta_actual, N_redes, "RESULTADOS_PRUEBA.txt");
    }
  
    /*
           evolucion_persona_a_persona(
        "ARCHIVOS_REDES/ER/ER_680.txt",
        "Resultados (PARTE 0)/Evolucion individual/ER/evoluciones/ER_680.txt",
        N_pasos,
        dt,
        K,
        betta,
        "Resultados (PARTE 0)/Evolucion individual/ER/velocidades/ER_680.txt"
    );
    
    */
    return 0;
}