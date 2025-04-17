#include "funciones.h"


int main() {
//   srand(time(NULL)); // Inicializa aleatoriedad
//   int N_sim=100;
//   int N_pasos=280;
//   double dt=0.1;
//   double betta=0.1;
//   int number_name=99;

//   muchas_simulaciones_ER(N_sim, N_pasos, dt,K, betta, number_name);


    srand(time(NULL)); // Inicializa aleatoriedad

    double K=10;         // a partir de ahora fijamos K=10
    int N_sim=200;
    int N_pasos=280;   
    double dt=0.1;
    double betta_max=1.0;
    double betta=0;
    double N, delta_betta;
    N=100;
    delta_betta=(betta_max -betta)/N;
    //empezamos con B=0 y realizamos los caluculos, aumentamos la betta un
    // delta_betta y reprtimos el proceso, y asi sucesivamente hasta B=1.

    //queremos ver cuál es la fracción de estados en los que obtenemos polarización.
    //Para ello, vamos a simular 200 veces para cada valor de betta para cada <k>.
    //las redes ER con <k>=6 son las que van de ER_200 a ER_399
    //las redes ER con <k>=14 son las que van de ER_400 a ER_599

    int number_name=399;
    int n_pol_6=0;
    int n_pol_14=0;
    

    //bucle para <k>=6
    //aqui vamos a hacer 200 simulaciones para cada valor de betta, y al final de cada simulacion
    //vamos a ver si la red es polarizada o no, y lo guardamos en un fichero.


    int=ult_num
    char* nom_carp_out = "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\polarizacion\\6";
    ult_num_out = obtener_siguiente_indice(nom_carp_out);
    char filename_out[512];
    sprintf(filename_out, "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\polarizacion\\6\\ER6_%d.txt", ult_num_out);
    //ahora tenemos un fichero donde vamos a guardar los resultados de la polarizacion de cada red.


    int num_polarizadas=0;
    char* nom_carp_in
    //bucle para <k>=6

    for(int i=0;i<N;i++) {
        muchas_simulaciones_ER(N_sim, N_pasos, dt,K, betta, number_name);
        nom_carp_in = "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\Evolucion temporal (promedio y desvest)\\ER";
        ult_num_in = obtener_siguiente_indice(nom_carp_in);
        //me hago las 200 simulaciones para esa betta
        for(int j=0; j<N_sim;j++){
            //miro mis 200 ficheros y cuento cuantas son polarizadas
            //paso el archivo a la funcion es_polarizada y me devuelve un 1 o un 0.
            
            char fich_in[512];
            sprintf("C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\Evolucion temporal (promedio y desvest)\\ER\\ER_%d.txt", ult_num_in-j);
            if (es_polarizada(fich_in)) {
                num_polarizadas++;
            }

        }

        //ahora guardo el resultado en el fichero de salida
        fopen(filename_out,"a");
        fprintf(filename_out, "%f\t %d\n", betta, num_polarizadas);
        fclose(filename_out);

        betta=betta+delta_betta;
    }
    //ahora tenemos un fichero donde tenemos la betta y el numero de redes polarizadas para cada betta.



    number_name=599;
    betta=0;
    nom_carp_out = "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\polarizacion\\14";
    ult_num_out = obtener_siguiente_indice(nom_carp_out);
    sprintf(filename_out, "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\polarizacion\\14\\ER14_%d.txt", ult_num_out);
    //ahora tenemos un fichero donde vamos a guardar los resultados de la polarizacion de cada red.


    num_polarizadas=0;
    //bucle para <k>=6

    for(int i=0;i<N;i++) {
        muchas_simulaciones_ER(N_sim, N_pasos, dt,K, betta, number_name);
        nom_carp_in = "C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\Evolucion temporal (promedio y desvest)\\ER";
        ult_num_in = obtener_siguiente_indice(nom_carp_in);
        //me hago las 200 simulaciones para esa betta
        for(int j=0; j<N_sim;j++){
            //miro mis 200 ficheros y cuento cuantas son polarizadas
            //paso el archivo a la funcion es_polarizada y me devuelve un 1 o un 0.
            
            char fich_in[512];
            sprintf("C:\\Users\\HP\\Desktop\\FISICA\\3 (2024-2025)\\segundo cuatri\\caos\\trabajo\\CaOtIcOs\\Resultados (PARTE 4)\\Evolucion temporal (promedio y desvest)\\ER\\ER_%d.txt", ult_num_in-j);
            if (es_polarizada(fich_in)) {
                num_polarizadas++;
            }

        }

        //ahora guardo el resultado en el fichero de salida
        fopen(filename_out,"a");
        fprintf(filename_out, "%f\t %d\n", betta, num_polarizadas);
        fclose(filename_out);

        betta=betta+delta_betta;
    }



    //bueno si estoy quisiera funcionar seria top, pero como no es el caso pues xd, contianuamos.

}