clear
// ====================================================================
// MODELO CINÉTICO DE DIFUSION + REACCION (Vectorizado)
// ====================================================================

// 1. Función para precalcular la matriz de coeficientes cinéticos K
function K = build_K(v, A, vo)
    K = zeros(v, v);
    for r = 1:v
        for c = 1:v
            term1 = 1/(r^(1/3));
            term2 = 1/(c^(1/3));
            term3 = r^(1/3);
            term4 = c^(1/3);
            K(r, c) =(term3 + term4)*( ( (A * (term1 + term2))^(-1) + (vo * (term3 + term4))^(-1) )^(-1));
        end
    end
endfunction

// 2. Función f2d Vectorizada (Primera etapa)
function [dif] = f2d_vec(y, K_mat, v, ycat, kred)
    dif = zeros(v, 1);
    sum2_vector = K_mat * y; // Vectorización del consumo
    dif(1)=-kred*ycat;
    for j = 1:v 
        sum1 = 0;
        for n = 1:j-1
            sum1 = sum1 + K_mat(j-n, n) * y(j-n) * y(n);
        end 
        if j==1 then
            sum1=kred*ycat;
            dif(j+1)=sum1- y(j) * sum2_vector(j);
        else
            dif(j+1) = 0.5 * sum1 - y(j) * sum2_vector(j);
        end
        
    end
endfunction

// 3. Función f3d Vectorizada (Segunda etapa, considera ultimoactivo)
function [dif] = f3d_vec(y, K_mat, v, ultimoactivo, ycat, kred)
    dif = zeros(v, 1);
    // Vectorización del consumo limitando hasta ultimoactivo
    sum2_vector = K_mat(1:v, 1:ultimoactivo) * y(1:ultimoactivo); 
    dif(1)=-kred*ycat;
    for j = 1:v 
        sum1 = 0;
        for n = 1:j-1
            sum1 = sum1 + K_mat(j-n, n) * y(j-n) * y(n);
        end 
        if j==1 then
            sum1=kred*ycat;
            dif(j+1)=sum1- y(j) * sum2_vector(j);
        else
            dif(j+1) = 0.5 * sum1 - y(j) * sum2_vector(j);
        end
    end
endfunction

// ====================================================================
// PARÁMETROS E INICIALIZACIÓN
// ====================================================================
t_inicio=getdate("s"); 
timer();
cpu_acumulado = 0;
v = 15;          // Tamaño de partícula mayor inicial
A = 1e-9;        // Constante del sistema
no = 500000;     // Condición inicial
h = 1;           // Paso de tiempo
t = 0;           // Tiempo inicial
tf = 1/(no*A);   // Tiempo final de la primera etapa
n = round(tf/h); // Iteraciones primera etapa
vo=0.2;
kred=0.01;
// Precalculamos la matriz K para la malla inicial
K_mat = build_K(v, A, vo);
// Inicialización de matrices de resultados
L = zeros(n, v+1);  
diferencial = zeros(n, v+1);
masatotal = zeros(1, n);
iterfinal=50000; //Tiempo final del desarrollo

ycationes(1) = no;  
L(1, 1) = t;
diferencial(1, 1) = t;
L(1, 2:v) = 0; 
contIter=0;       
tiempoAnt=0;

// ====================================================================
// ETAPA 1: Método Runge-Kutta 4 (Malla Fija)
// ====================================================================
disp("Iniciando Simulacion...");

// ====================================================================
// ETAPA 2: Método Runge-Kutta 4 (Malla Dinámica)
// ====================================================================
tol = 10^(-5); // Tolerancia para la actividad

for i = 2:iterfinal
    // Vectorización de la actividad
    actividad = L(i-1, 2:v+1) >= tol;
    
    // Encontrar ultimoactivo
    idx_falso = find(~actividad, 1);
    if isempty(idx_falso) then
        ultimoactivo = length(actividad);
    else
        ultimoactivo = idx_falso - 1;
        if ultimoactivo == 0 then ultimoactivo = length(actividad); end
    end
    
    tamañomax = ultimoactivo * 2;
    verificacion = %T;
    
    // Expansión dinámica de la malla
    while verificacion == %T 
        if tamañomax > v then
            v = v + 1;
            L(:, v+1) = 0;
            diferencial(:, v+1) = 0;
            // Al agrandar v, reconstruimos la matriz K_mat
            K_mat = build_K(v, A); 
            verificacion = %T;
        else
            verificacion = %F;
        end
    end
    
    y_n = L(i-1, 2:v+1)'; 
    ycat=ycationes(i-1);
    if ycat<0 then
        ycat=0;
    end
    // RK4 para Etapa 2
    k1 = f3d_vec(y_n, K_mat, v, ultimoactivo, ycat, kred);
    y_naux(1:v)=y_n(:) + k1(2:v+1)*h/2;
    ycataux=ycat+k1(1)*h/2;
    k2 = f3d_vec(y_naux, K_mat, v, ultimoactivo, ycataux, kred);
    y_naux(:)=y_n(:) + k2(2:v+1)*h/2;
    ycataux=ycat+k2(1)*h/2;
    k3 = f3d_vec(y_naux, K_mat, v, ultimoactivo, ycataux, kred);
    y_naux(:)=y_n(:) + k3(2:v+1)*h;
    ycataux=ycat+k3(1)*h;
    k4 = f3d_vec(y_naux, K_mat, v, ultimoactivo, ycataux, kred);
    
    // Combinación
    L(i, 2:v+1) = y_n(1:v) + (h/6) * (k1(2:v+1) + 2*k2(2:v+1) + 2*k3(2:v+1) + k4(2:v+1));
    diferencial(i,2)=(k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6;
    diferencial(i, 3:v+2) = (k1(2:v+1) + 2*k2(2:v+1) + 2*k3(2:v+1) + k4(2:v+1)) / 6;
    ycationes(i)=ycat+h* (k1(1) + 2*k2(1) + 2*k3(1) + k4(1))/6;
    // Actualización tiempo y masa
    t = t + h;
    L(i, 1) = t;
    diferencial(i, 1) = t;
    masatotal(i) = sum(L(i, 2:v+1) .* (1:v)) + ycationes(i);
  
    // Monitorización
    if modulo(i, 1000) == 0 then
        mprintf("Iter %d - t=%.2f - Masa=%.2f - y1=%.4f\n", i, t, masatotal(i), L(i, 2));
        contIter=contIter+1;
        tiempoActual = getdate("s") - t_inicio; 
        cpu_delta = timer();
        cpu_acumulado = cpu_acumulado + cpu_delta;
        tiempoIter(contIter,1) = i;
        tiempoIter(contIter,2) = tiempoActual;
        tiempoIter(contIter,3) = tiempoActual - tiempoAnt;
        tiempoIter(contIter,4) = cpu_acumulado;
        tiempoIter(contIter,5) = cpu_delta;
        // Actualizamos para la próxima vuelta
        tiempoAnt = tiempoActual;
    end   
end

t_etapa1 = getdate("s") - t_inicio;
disp("Simualacion finalizada en " + string(t_etapa1) + " segundos.");
nombreCarpeta = "Lineal de "+string(no)+" partículas iniciales y de "+string(iterfinal)+" iteraciones don Difusion, R. Sup. y Reduccion";
ruta_base = "/home/bianca/Documentos/2026/Marzo/";
ruta_carpeta = ruta_base + nombreCarpeta;
mkdir(ruta_carpeta);
ResultadoTotal(:,1)=L(:,1);
ResultadoTotal(1:iterfinal,2:v+1)=L(1:iterfinal,2:v+1)*diag(1:v);
csvWrite(ResultadoTotal,ruta_carpeta+"/ResultadoTotal.csv")
ResultadoTotalConcentracion=L;
csvWrite(ResultadoTotalConcentracion,ruta_carpeta+"/ResultadoConcentracion.csv")
ResultadoDiferenciales=diferencial;
csvWrite(diferencial,ruta_carpeta+"/ResultadoDiferenciales.csv")
csvWrite(tiempoIter, ruta_carpeta+"/TiempoIteracion.csv")
csvWrite(ycationes, ruta_carpeta+"/ResultadoCatalizador.csv")
