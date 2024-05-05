#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define MAX_X 1000 // Valor máximo para la coordenada x
#define MAX_Y 1000 // Valor máximo para la coordenada y

// Estructura para representar un punto en el plano
typedef struct {
    int x, y;
} Punto;

//Funcion para obtener la regresion cuadrática
void regresionPolinomial(double x[], double y[], unsigned long n, double *a, double *b, double *c){
	const int grado = 2; // Cambiar el grado para cambiar el grado del polinomio resultante

    double matriz[grado + 1][grado + 1];
    double soluciones[grado + 1];

    double sumX = 0;
    double sumXY = 0;


    // Inicializamos el sistema de ecuaciones
    for (int i = 0; i < grado + 1; ++i) {
        sumXY = 0;
        for (int k = 0; k < n; ++k) {
            sumXY += pow(x[k], i) * y[k];
        }
        soluciones[i] = sumXY;

        for (int j = 0; j < grado + 1; ++j) {
            if (i == 0 && j == 0) {
                matriz[i][j] = n;
            } else {
                sumX = 0;
                for (int k = 0; k < n; ++k) {
                    sumX += pow(x[k], (j + i));
                }
                matriz[i][j] = sumX;
            }
        }
    }

    // Resolvemos el sistema de ecuaciones con el método de Gauss-Seidel
    int size = grado + 1;
    double EPSILON = 0.0000001;

    double S[size];
    double OS[size];
    double E[size];
    for (int i = 0; i < size; ++i) {
        S[i] = 0;
        OS[i] = 0;
        E[i] = 0;
    }

    int errorCheck = 0;
    int counter = 0;

    while (!errorCheck) {
        for (int i = 0; i < size; ++i) {
            double sum = 0;
            for (int j = 0; j < size; ++j) {
                if (j != i) {
                    sum += matriz[i][j] * S[j];
                }
            }
            OS[i] = S[i];
            S[i] = (soluciones[i] - sum) / matriz[i][i];
            E[i] = (S[i] - OS[i]);
        }

        errorCheck = 1;
        for (int i = 0; i < size; ++i) {
            if (fabs(E[i]) > EPSILON) {
                errorCheck = 0;
                break;
            }
        }
        counter++;
    }

    printf("------------------POLINOMIO---------------\n");
    for (int i = size - 1; i >= 1; --i) {
        printf("%fX^%d + ", S[i], i);
    }
    printf("%f\n", S[0]);
    
    *a=S[2];
    *b=S[1];
    *c=S[0];
    
}


void regresionLineal(double x[], double y[], int n, double *a, double *b){
	// Declaración de variables
    int i;
    double m, c, sumax, sumay, sumaxy, sumax2;

    // Calculamos las sumatorias
    sumax = sumay = sumaxy = sumax2 = 0;
    for (i = 0; i < n; i++) {
        sumaxy += x[i] * y[i];
        sumax2 += x[i] * x[i];
        sumax += x[i];
        sumay += y[i];
    }
    
    // Calculamos la pendiente (m) y la intersección (c)
    m = (n * sumaxy - sumax * sumay) / (n * sumax2 - sumax * sumax);
    c = (sumay - m * sumax) / n;
    
    *a = m;
    *b = c;

    // Mostramos la ecuación de la regresión lineal
    printf("Polinomio lineal: Y = %.2f * X + (%.2f)\n\n", m, c);
}


// Función de comparación para ordenar los puntos por coordenada x
int comparar_x(const void *a, const void *b) {
    Punto *p1 = (Punto *)a;
    Punto *p2 = (Punto *)b;
    return p1->x - p2->x;
}

// Función de comparación para ordenar los puntos por coordenada y
int comparar_y(const void *a, const void *b) {
    Punto *p1 = (Punto *)a;
    Punto *p2 = (Punto *)b;
    return p1->y - p2->y;
}

// Función para calcular la distancia entre dos puntos
double distancia(Punto p1, Punto p2) {
    return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
}

// Función para fusionar dos arreglos ordenados
void fusionar(Punto arr[], int izquierda, int medio, int derecha, int(*comparar)(const void *, const void *)) {
    int tamano_izquierda = medio - izquierda + 1;
    int tamano_derecha = derecha - medio;

    Punto *arreglo_izquierda = malloc(tamano_izquierda * sizeof(Punto));
    Punto *arreglo_derecha = malloc(tamano_derecha * sizeof(Punto));

    // Copiar los datos en los subarreglos
    for (int i = 0; i < tamano_izquierda; i++) {
        arreglo_izquierda[i] = arr[izquierda + i];
    }
    for (int i = 0; i < tamano_derecha; i++) {
        arreglo_derecha[i] = arr[medio + 1 + i];
    }

    int i = 0, j = 0, k = izquierda;
    while (i < tamano_izquierda && j < tamano_derecha) {
        if (comparar(&arreglo_izquierda[i], &arreglo_derecha[j]) <= 0) {
            arr[k++] = arreglo_izquierda[i++];
        } else {
            arr[k++] = arreglo_derecha[j++];
        }
    }

    while (i < tamano_izquierda) {
        arr[k++] = arreglo_izquierda[i++];
    }
    while (j < tamano_derecha) {
        arr[k++] = arreglo_derecha[j++];
    }

    free(arreglo_izquierda);
    free(arreglo_derecha);
}

// Función para ordenar un arreglo de puntos por coordenada y
void merge_sort(Punto arr[], int izquierda, int derecha, int(*comparar)(const void *, const void *)) {
    if (izquierda < derecha) {
        int medio = izquierda + (derecha - izquierda) / 2;

        merge_sort(arr, izquierda, medio, comparar);
        merge_sort(arr, medio + 1, derecha, comparar);

        fusionar(arr, izquierda, medio, derecha, comparar);
    }
}

// Función para calcular la distancia mínima entre pares de puntos utilizando el algoritmo de divide y vencerás
double parMasCercano(Punto puntos[], int n, unsigned long long *tiempo_total) {
    if (n <= 1) {
        return INFINITY; // No hay suficientes puntos para calcular la distancia
    }

    struct timespec inicio, fin;
    clock_gettime(CLOCK_MONOTONIC, &inicio);

    int medio = n / 2;
    Punto puntoMedio = puntos[medio];

    // Dividir y conquistar
    double dl = parMasCercano(puntos, medio, tiempo_total);
    double dr = parMasCercano(puntos + medio, n - medio, tiempo_total);

    double d = fmin(dl, dr);

    // Crear un arreglo para almacenar puntos dentro del área de exploración
    Punto *tira = malloc(n * sizeof(Punto));
    int j = 0;
    for (int i = 0; i < n; i++) {
        if (abs(puntos[i].x - puntoMedio.x) < d) {
            tira[j] = puntos[i];
            j++;
        }
    }

    // Ordenar la tira de puntos por coordenada y
    merge_sort(tira, 0, j - 1, comparar_y);

    // Verificar puntos dentro de los 7 vecinos más cercanos en la tira
    for (int i = 0; i < j; i++) {
        for (int k = i + 1; k < j && (tira[k].y - tira[i].y) < d; k++) {
            double dist = distancia(tira[i], tira[k]);
            if (dist < d) {
                d = dist;
            }
        }
    }

    free(tira);

    clock_gettime(CLOCK_MONOTONIC, &fin);

    *tiempo_total += (fin.tv_sec - inicio.tv_sec) * 1000000000 + (fin.tv_nsec - inicio.tv_nsec);

    return d;
}

void registrarResultado(int n, unsigned long long tiempo) {
    FILE *archivo = fopen("datosfuerzabruta.txt", "a");
    if (archivo == NULL) {
        printf("Error al abrir el archivo.\n");
        return;
    }
    fprintf(archivo, "%d %llu\n", n, tiempo);
    fclose(archivo);
}

int main(void) {
	
	double a=0,b=0,c=0; // Variables para los coeficientes de la regresión cuadrática
    double a1=0,b1=0,a2=0,b2=0; // Variables para los coeficientes de la regresión lineal
    
	FILE *archivo; 
    archivo = fopen("datosfuerza.txt", "r");
	
    // Verificar si el archivo se abrió correctamente
    if (archivo == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }
    
    int max_size = 10000; // Cambia esto al tamaño máximo que esperas para los datos
    double puntos[max_size];
    double tiempos[max_size];
    int index = 0;

    // Leer los valores del archivo
    while (fscanf(archivo, "%lf %lf", &puntos[index], &tiempos[index]) == 2) {
        index++;
        // Verificar si hemos alcanzado el límite máximo de almacenamiento
        if (index >= max_size) {
            printf("¡Se alcanzó el límite máximo de almacenamiento!\n");
            break;
        }
    }

    // Cerrar el archivo
    fclose(archivo);
    
    
    regresionPolinomial(puntos,tiempos,index,&a,&b,&c); //regresion cuadrática
    
    
    FILE *archivo1; 
    archivo1 = fopen("datosdivide.txt", "r");

    // Verificar si el archivo se abrió correctamente
    if (archivo1 == NULL) {
        perror("Error al abrir el archivo");
        return 1;
    }

    double puntos1[max_size];
    double tiempos1[max_size];
    int index1 = 0;

    // Leer los valores del archivo
    while (fscanf(archivo, "%lf %lf", &puntos1[index1], &tiempos1[index1]) == 2) {
        index1++;
        // Verificar si hemos alcanzado el límite máximo de almacenamiento
        if (index1 >= max_size) {
            printf("¡Se alcanzó el límite máximo de almacenamiento!\n");
            break;
        }
    }

    // Cerrar el archivo
    fclose(archivo);
    
    
    // Realizar la regresión lineal
    regresionLineal(puntos1, tiempos1, index1, &a1, &b1);
	

    /*
    int umbral = 89; // Umbral para decidir qué algoritmo usar
    
    unsigned long long tiempoDivide = 0;

	for (int n = 3; n <= 200; n += 1) {
        // Repetir el proceso cinco veces para cada tamaño de punto
        for (int rep = 0; rep < 3; rep++) {
        	tiempoDivide=0;
        	
    		// Inicializar puntos en un arreglo
		    Punto *puntos = malloc(n * sizeof(Punto));
		    for (int i = 0; i < n; i++) {
		        puntos[i].x = rand() % (MAX_X + 1);
		        puntos[i].y = rand() % (MAX_Y + 1);
		    }
		
		    double distanciaMinima=0;
		
		    // Decidir qué algoritmo utilizar
		    if (n <= umbral) {
		        printf("Utilizando fuerza bruta para %d puntos.\n", n);
		        struct timespec inicio, fin;
	        	// Inicia la medición del tiempo
		        clock_gettime(CLOCK_MONOTONIC, &inicio);
		
		        double dm = INFINITY; // Inicializar dm con un valor infinito
	            double dt;
	            int it;
	            int jt;
	            for (int i = 0; i < n; i++) {
	                for (int j = 0; j < n; j++) {
	                    if (i != j) { // Evitar comparar un punto consigo mismo
	                        dt = sqrt(pow(puntos[j].x - puntos[i].x, 2) + pow(puntos[j].y - puntos[i].y, 2));
	                        if (dt < dm) {
	                            dm = dt;
	                            it = i;
	                            jt = j;
	                        }
	                    }
	                }
	            }
		        // Finaliza la medición del tiempo
		        clock_gettime(CLOCK_MONOTONIC, &fin);
		        unsigned long long tiempoFuerza = (fin.tv_sec - inicio.tv_sec) * 1000000000 + fin.tv_nsec - inicio.tv_nsec;
		        printf("Tiempo medido: %llu nanosegundos\n", tiempoFuerza);
		        registrarResultado(n, tiempoFuerza);   
		    } else {
		        printf("Utilizando divide y vencerás para %d puntos.\n", n);
		        merge_sort(puntos, 0, n - 1, comparar_x);
		        distanciaMinima = parMasCercano(puntos, n, &tiempoDivide);
		        free(puntos);
		        registrarResultado(n, tiempoDivide);
		        printf("Tiempo medido: %llu nanosegundos\n", tiempoDivide);
		    }
		}
	}
	
	*/
	
    
    FILE *script = fopen("grafico.gp", "w");
	if (script == NULL) {
	    perror("Error al abrir el archivo de script de Gnuplot");
	    return 1;
	}
	
	
	// Escritura del contenido del script
	fprintf(script, "# Establecer el título y etiquetas de los ejes\n");
	fprintf(script, "set title 'Gráfico con lumbral establecido'\n");
	fprintf(script, "set xlabel 'Cantidad de puntos'\n");
	fprintf(script, "set ylabel 'Tiempo (nanosegundos)'\n");
	fprintf(script, "# Graficar los puntos desde el archivo datosfuerza.txt, datosdivide.txt y datos.txt, y ajustar las funciones\n");
	fprintf(script, "plot 'datosfuerza.txt' with points pointtype 7 pointsize 0.7 linecolor rgb 'light-green' title 'Puntos Fuerza Bruta', \
	    'datosdivide.txt' with points pointtype 7 pointsize 0.7 linecolor rgb 'light-blue' title 'Puntos Divide y vencerás', \
	    %f*x**2 + %f*x + %f title 'Regresión Cuadrática' linewidth 1 linecolor 'green', \
	    %f*x + %f title 'Regresión Lineal' linewidth 1 linecolor 'blue', \
	    'datoslumbral.txt' with points pointtype 7 pointsize 0.3 linecolor rgb 'light-red' title 'Ajuste de lumbral' ", a, b, c, a1, b1);

	// Cerrar el archivo de script
	fclose(script);
	
	// Ejecutar Gnuplot con el script
	FILE *gnuplot_pipe = popen("gnuplot -persist", "w");
	if (gnuplot_pipe == NULL) {
	    perror("Error al abrir el pipe de Gnuplot");
	    return 1;
	}
	// Cargar el script de Gnuplot
	fprintf(gnuplot_pipe, "load 'grafico.gp'\n");
	// Cerrar el pipe de Gnuplot
	pclose(gnuplot_pipe);
    
	
    return 0;
}

