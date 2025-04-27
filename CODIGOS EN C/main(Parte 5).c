#include "funciones.h"
/*  Notad que, si el grado es lo suficientemente pequeño (y por tanto el horizonte
    de información), puede aparecer polarización incluso en ausencia de homofilia
    (beta = 0). Estudiad sistemáticamente esta aparición con redes ER de grado bajo,
    y con redes WS de grado k= 10 con diversas probabilidades de rewiring
    (reproduciendo cualitativamente la figura 4 del paper, pero con fracciones de
    polarización).*/    


    #include "funciones.h"


int main() {
        //hola caracola
        srand(time(NULL)); // Inicializa aleatoriedad
        double dt=0.1;
        int N_pasos=50;
        double K=10;
        int N_redes=100;

        int rede_ini=0;
    
        double betta_= 0 ;
        double k_inicial=2;
        double k_final=8;
        double delta_k=0.1;

        double k_actual=k_inicial;

        while(k_actual<k_final){
            frac_polarizado_apartado5(N_redes, rede_ini,k_actual, N_pasos,dt);
            calcular_fraccion_polarizados_apartado5(k_actual, N_redes, "Resultados.txt");
            k_actual=k_actual+delta_k;
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
    