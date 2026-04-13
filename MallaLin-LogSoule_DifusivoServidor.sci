clear

function K = build_K(v, A)
    K = zeros(v, v);
    for r = 1:v
        for c = 1:v
            term1 = 1/(r^(1/3));
            term2 = 1/(c^(1/3));
            term3 = r^(1/3);
            term4 = c^(1/3);
            K(r, c) = A * (term1 + term2) * (term3 + term4);
        end
    end
endfunction

function K = build_KLog(nmalla, A, valMedio)
    K = zeros(nmalla, nmalla);
    for r = 1:nmalla
        for c = 1:nmalla
            //pause
            term1 = 1/(valMedio(c)^(1/3));
            term2 = 1/(valMedio(r)^(1/3));
            term3 = valMedio(c)^(1/3);
            term4 = valMedio(r)^(1/3);
            K(r, c) = A * (term1 + term2) * (term3 + term4);
        end
    end
endfunction

function [dif] = f2d(i, L, v,Coef, ultimoactivo) //funcion para el inicio del programa
    dif = zeros(1, v);
    y(:)=L(i-1,2:ultimoactivo+1);
    // Vectorización del consumo limitando hasta ultimoactivo
    sum2_vector = Coef(1:v, 1:ultimoactivo) * y(:); 
    for j = 1:v 
        sum1 = 0;
        for n = 1:j-1
            sum1 = sum1 + Coef(j-n, n) * L(i-1,j-n+1) * L(i-1,n+1);
        end 
        dif(1,j) = 0.5 * sum1 - L(i-1,j+1) * sum2_vector(j);
    end
endfunction



//Derivada para la parte logaritmica
function [diferencialLog] = funcionLogaritmica(i, L2, ultimoactivo, linealizacion, valMedio, diflimites, nmalla, limites, Coef) 
    // Vector de salida  
    q=zeros(nmalla,1);
    b=zeros(nmalla,1);
    for l=linealizacion+1:nmalla
        if l==nmalla || L2(i-1,l+2)==0 then
            numerador_q = log10(L2(i-1,l+1)/(limites(l,2)-limites(l-1,2))) - log10(L2(i-1,l-1)/(limites(l-1,1)-limites(l-2,1)));
            denominador_q=(log10(valMedio(l))-log10(valMedio(l-2)));
            q(l,1)=numerador_q/denominador_q;
        else
        numerador_q = log10(L2(i-1,l+2)/((limites(l+1,2))-limites(l,2))) - log10(L2(i-1,l)/(limites(l,1)-limites(l-1,1)));
        denominador_q=(log10(valMedio(l+1))-log10(valMedio(l-1)));
        q(l,1)=numerador_q/denominador_q;    
        end
        b(l,1) = (q(l,1)+1)*L2(i-1,l+1)*(limites(l,1)^q(l,1))/((limites(l,2)+1)^(q(l,1)+1)-limites(l,1)^(q(l,1)+1));
    end
    for j = 1:nmalla
        sum1 = 0;
        sum2 = 0;
        sumL = 0;
        sumLS = 0;
        sumLI = 0;
        sumMM = 0;
        sumMMLI=0;
        sumMMLS=0;
        sumML=0;
        sum2vector=0;
        if j>=linealizacion+1 then
            
            if limites(j,1)<(linealizacion*2) then
                if limites(j,2)>linealizacion*2 then
                    val=linealizacion*2; 
                else
                    val=limites(j,2);
                end
                
            // Primer sumatorio - versión completamente corregida
                for l = limites(j,1):val
                    for n = l-linealizacion:linealizacion
                         sumL = sumL + Coef(l-n,n) * L2(i-1, l-n+1) * L2(i-1, n+1);
                        
                    end
                end
            end
            
                if j > ultimoactivo then
                    valor = ultimoactivo;
                else
                    valor = j;
                end
                
                for l = linealizacion+1:valor
                    for k = 1:linealizacion
                        tamañoInf = limites(l,1) + valMedio(k);
                        tamañoSup = limites(l,2)+1 + valMedio(k);
                        if tamañoInf >= limites(j,1) & tamañoSup < limites(j,2)+1 then
                            sumML = sumML +Coef(l,k) * L2(i-1, k+1) * L2(i-1, l+1);
                            
                        elseif tamañoInf < limites(j,2)+1 & tamañoSup >= limites(j,2)+1 then
                            valInf = limites(l,1);
                            valSup = limites(j,2)+1-valMedio(k);
                            if q(l,1)<> -1 then
                                Ninf = b(l,1)*limites(l,1)*((valSup/limites(l,1))^(q(l,1)+1)-(valInf/limites(l,1))^(q(l,1)+1))/(q(l,1)+1);
                            else
                                Ninf = b(l,1)*limites(l,1)*(log(valSup/limites(l,1))-log(valInf/limites(l,1)));
                            end
                            sumLI = sumLI +  Coef(l,k)*Ninf*L2(i-1,k+1);
                           
                        elseif tamañoInf <= limites(j,1) & tamañoSup >= limites(j,1) then
                            valInf = limites(j-1,2)+1-valMedio(k);
                            valSup = limites(l,2)+1;
                            if q(l,1)<> -1 then
                                Nsup = b(l,1)*limites(l,1)*((valSup/limites(l,1))^(q(l,1)+1)-(valInf/limites(l,1))^(q(l,1)+1))/(q(l,1)+1);
                            else
                                Nsup = b(l,1)*limites(l,1)*(log(valSup/limites(l,1))-log(valInf/limites(l,1)));
                            end
                            sumLS = sumLS + Coef(l,k)*Nsup*L2(i-1,k+1);
                        end
                    end
                
                 end
                 for l = linealizacion+1:valor   
                    for k = l:valor
                        if k == l then
                            menor = k;
                            mayor = l;
                        else
                            menor = min(l,k);
                            mayor = max(l,k);
                        end
                        
                        tamañoInf = limites(mayor,1) + valMedio(menor);
                        tamañoSup = limites(mayor,2)+1 + valMedio(menor);
                        
                        if tamañoInf >= limites(j,1) & tamañoSup < limites(j,2)+1 then
                            sumMM = sumMM + Coef(menor,mayor) * L2(i-1, menor+1) * L2(i-1, mayor+1);
                            
                        elseif tamañoInf < limites(j,2)+1 & tamañoSup >= limites(j,2)+1 then
                            valInf = limites(mayor,1);
                            valSup = limites(j,2)+1-valMedio(menor);
                            if q(mayor,1)<> -1 then
                                Ninf = b(mayor,1)*limites(mayor,1)*((valSup/limites(mayor,1))^(q(mayor,1)+1)-(valInf/limites(mayor,1))^(q(mayor,1)+1))/(q(mayor,1)+1);
                            else
                                Ninf = b(mayor,1)*limites(mayor,1)*(log(valSup/limites(mayor,1))-log(valInf/limites(mayor,1)));
                            end
                            sumMMLI = sumMMLI +Coef(menor,mayor)*Ninf*L2(i-1,menor+1);
                            
                        elseif tamañoInf <= limites(j,1) & tamañoSup >= limites(j,1) then
                            valInf = limites(j-1,2)+1-valMedio(menor);
                            valSup = limites(mayor,2)+1;
                            
                            if q(l,1)<> -1 then
                                Nsup = b(mayor,1)*limites(mayor,1)*((valSup/limites(mayor,1))^(q(mayor,1)+1)-(valInf/limites(mayor,1))^(q(mayor,1)+1))/(q(mayor,1)+1);
                            else
                                Nsup = b(mayor,1)*limites(mayor,1)*(log(valSup/limites(mayor,1))-log(valInf/limites(mayor,1)));
                            end
                            sumMMLS = sumMMLS + Coef(menor,mayor)*Nsup*L2(i-1,menor+1);
                            
                        end
                   
                end
            end
            // Segundo sumatorio - versión completamente corregida
            if j<= ultimoactivo then
                for n = 1:ultimoactivo
                   sum2vector = Coef(j,1:ultimoactivo) .* L2(i-1, 2:ultimoactivo+1); 
                    sum2=sum(sum2vector);
                end
            end
            
            sum1 = 0.5 * sumL + sumLS + sumLI + sumMM + sumMMLI + sumMMLS + 0.5*sumML;
            diferencialLog(1,j) = (sum1 - L2(i-1, j+1) * sum2);
            diferencialLog(2,j)= sumL*0.5;
            diferencialLog(3,j)=sumLS;
            diferencialLog(4,j)=sumLI;
            diferencialLog(5,j)=sumMM;
            diferencialLog(6,j)=sumMMLS;
            diferencialLog(7,j)=sumMMLI;
            diferencialLog(8,j)=sumML*0.5;
            diferencialLog(9,j)=-sum2 *L2(i-1, j+1);
       else
            // Primer sumatorio - versión completamente corregida
            if j<>1 then 
                for n = 1:j-1
                    sum1 = sum1 + Coef(j-n,n) * L2(i-1, j-n+1) * L2(i-1, n+1);
                    sumL=sum1;
                end 
            end
            
            // Segundo sumatorio - versión completamente corregida
            sum2vector = Coef(j,1:ultimoactivo) .* L2(i-1, 2:ultimoactivo+1); 
            sum2=sum(sum2vector);
            
            diferencialLog(1,j) = (0.5 * sum1 - L2(i-1, j+1) * sum2);
            diferencialLog(2,j)= 0.5*sumL;
            diferencialLog(3,j)=sumLS;
            diferencialLog(4,j)=sumLI;
            diferencialLog(5,j)=sumMM;
            diferencialLog(6,j)=sumMMLS;
            diferencialLog(7,j)=sumMMLI;
            diferencialLog(8,j)=sumML;
            diferencialLog(9,j)=-sum2 *L2(i-1, j+1);
            
       end
    end
    
endfunction

//funcion que permite la distribución de la primera parte
function [normalizacion] = fNorm(ci, v, limites, nmalla, linealizacion)
    sum = 0;
    cont = 0;
    for i = 1:nmalla
        if i<=linealizacion then
            normalizacion(i) = ci(i);
        else
            num2= limites(i,2);
            num1=limites(i,1);
            for j=num1:num2
                sum= ci(j)+sum;
            end            
            normalizacion(i)=sum;
        end
        sum = 0;
    end
endfunction

//Definicion de actividad
function [actividad]=act(i,L2, tol,nmalla)    
    for j=1:nmalla
        if  L2(i,j+1)>= tol then
            actividad(j)=%T;
        else 
            actividad(j)=%F;
        end 
    end
endfunction

// --- Configuración inicial ---
t_inicio= getdate("s");
timer();
cpu_acumulado=0;  
v = 15;       //Tamaño de particula mayor
A = 1e-9;       // Constante del sistema
no = 500000;       // Condición inicial, cantidad de particulas unitariasí totales
h = 1;       // Paso de tiempo
t = 0;          // Tiempo inicial
tf = 1/(no*A);  // Tiempo final, valor recopilado de paper que es hasta donde el sistema evita la perdida de masa, a partir de ahí hay que forzarlo a no perder masa.
n = 5000;// Iteraciones
tol=10^(-5);
// --- Inicialización ---
L = zeros(n, v+1);  
L(1, 1) = t;
diferencial(1,1)=t;
L(1, 2) = no;  
for i = 2:v
    L(1, i+1) = 0; 
end           
ultimoactivo=v;
Coef=build_K(v,A);

contIter=0;
tiempoAnt=0;
// --- Método Runge-Kutta ---
for i = 2:n+1
    actividad = L(i-1, 2:v+1) >= tol;
    
    // Encontrar ultimoactivo
    idx_falso = find(~actividad, 1);
    if isempty(idx_falso) then
        ultimoactivo = length(actividad);
    else
        ultimoactivo = idx_falso - 1;
        if ultimoactivo == 0 || ultimoactivo ==1 then ultimoactivo = length(actividad); end
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
            Coef = build_K(v, A); 
            verificacion = %T;
        else
            verificacion = %F;
        end
    end
    // Paso 1
    k1 = f2d(i, L, v,Coef, ultimoactivo);
    Laux(i-1, 2:v+1) = L(i-1, 2:v+1) + k1(1,1:v) * h/2;
    
    // Paso 2
    k2 = f2d(i, Laux, v,Coef, ultimoactivo);
    //Laux = L;
    Laux(i-1, 2:v+1) = L(i-1, 2:v+1) + k2(1,1:v) * h/2;
    // Paso 3
    k3 = f2d(i, Laux, v,Coef, ultimoactivo);
    //Laux = L;
    Laux(i-1, 2:v+1) = L(i-1, 2:v+1) + k3(1,1:v) * h;
    // Paso 4
    k4 = f2d(i, Laux, v,Coef, ultimoactivo);
    
    // Combinación
    L(i, 2:v+1) = L(i-1, 2:v+1) + h/6 * (k1(1,1:v) + 2*k2(1,1:v) + 2*k3(1,1:v) + k4(1,1:v));
    diferencial(i-1,2:v+1)= h*(k1(1,1:v) + 2*k2(1,1:v) + 2*k3(1,1:v) + k4(1,1:v))/6;
    // Actualización tiempo
    t = t + h;
    L(i, 1) = t;
    diferencial(i-1,1)=t
    masatotal(i) = sum(L(i, 2:v+1).*(1:v));
  
    // Monitorización
    if modulo(i, 1000) == 0 then
        total = sum(L(i, 2:v+1).*(1:v));
        mprintf("Iter %d - t=%.2f - Masa=%.2f - y1=%.4f\n", i, t, total, L(i, 2));
        contIter=contIter+1;
        tiempoActual = getdate("s") - t_inicio; 
        cpu_delta = timer(); // <-- Tiempo CPU de estas últimas 1000 iteraciones
        cpu_acumulado = cpu_acumulado + cpu_delta; // <-- Tiempo CPU total acumulado
        tiempoIter(contIter,1) = i;
        tiempoIter(contIter,2) = tiempoActual;
        tiempoIter(contIter,3) = tiempoActual - tiempoAnt;
        tiempoIter(contIter,4) = cpu_acumulado; 
        tiempoIter(contIter,5) = cpu_delta;
       
        // Actualizamos para la próxima vuelta
        tiempoAnt = tiempoActual;
    end   
end



// Asignar valores de L a ci
for i = 1:v
    ci(i) = L(min(20000,n), i+1);
end

h=1;
z=1.2; // base logaritmica del espaciado
//linealizacion=v
iterfinal=50000;
linealizacion=min(find(ci<(max(ci)/10000)))+1; //Define donde comienza la parte logartimica
var =0;
limites(:,1)=1:linealizacion; //Matriz que contiene los valores de los tamanios de limites de cada malla.
limites(:,2)=0;
cont= linealizacion;
limites(cont +1 ,1)= linealizacion + 1 ;
while var==0
    cont = cont + 1;
    limites(cont,2)=round(limites(cont,1) + (z^(cont-linealizacion)))-1;
    if limites(cont,2)>v then
        var =1;
        limites(cont,2)=v
    else
        limites(cont+1,1)= limites(cont,2)+1;
    end
end
nmalla=length(limites)/2;
ciN = fNorm(ci, v, limites, nmalla, linealizacion); // normaliza las condiciones iniciales
tol=10^(-5); //tolerancia para la actividad
limites(nmalla,2)=round(limites(cont,1) + (z^(cont-linealizacion)))-1;
for i=1:nmalla
    if i<=linealizacion
        valMedio(i)=limites(i,1);
        difLimites(i)=limites(i,1);
    else 
        valMedio(i)=(limites(i,1)+limites(i,2))/2;
        difLimites(i)=limites(i,2)+1-limites(i,1);
    end
end 
Coef=build_KLog(nmalla, A, valMedio);
verificacion=%T
L2(1,1)=tf;
for i=2:length(ciN)+1
    L2(1,i)=ciN(i-1);
end
i=1;
for m=1:nmalla
   total = total +L2(i, m+1)*valMedio(m);
end
while i<>iterfinal-n
    actividad = act(i,L2, tol,nmalla);
    ultimoactivo= find(actividad==%F,1)-1;
    if isempty(ultimoactivo)|| ultimoactivo==0 then
        ultimoactivo=length(actividad);
        limites(ultimoactivo+1,1)=limites(ultimoactivo,2)+1;
        limites(ultimoactivo+1,2)=round(limites(ultimoactivo+1,1) + (z^(nmalla+1-linealizacion)))-1;
        ciN(nmalla+1)=0;
        nmalla=nmalla+1;
        L2(:,nmalla+1)=0;
        verificacion=%T;
        valMedio(nmalla)=(limites(nmalla,1)+limites(nmalla,2))/2;
        difLimites(nmalla)=limites(nmalla,2)+1-limites(nmalla,1);
        Coef=build_KLog(nmalla, A, valMedio);
    end
    tamañomax= limites(ultimoactivo+1,2) + (limites(ultimoactivo+1,2)+limites(ultimoactivo+1,1))/2;
    verificacion=%T;
    while verificacion==%T 
        if tamañomax>limites(nmalla,2) then
            limites(nmalla+1,1)=limites(nmalla,2)+1;
            limites(nmalla+1,2)=round(limites(cont,1) + (z^(nmalla+1-linealizacion)))-1;
            ciN(nmalla+1)=0;
            nmalla=nmalla+1;
            L2(:,nmalla+1)=0;
            verificacion=%T;
            valMedio(nmalla)=(limites(nmalla,1)+limites(nmalla,2))/2;
            difLimites(nmalla)=limites(nmalla,2)+1-limites(nmalla,1);
            Coef=build_KLog(nmalla, A, valMedio);
        else
            verificacion=%F;
        end
    end
    i=i+1;
    L2(i,1)=t;
    diferencial(i+1999,1)=t;
    diferencialLineal(i,1)=t;
    diferencialLS(i,1)=t;
    diferencialLI(i,1)=t;
    diferencialMM(i,1)=t;
    diferencialMMLS(i,1)=t;
    diferencialMMLI(i,1)=t;
    diferencialML(i,1)=t;
    diferencialSalida(i,1)=t;
    t=t+h;
    // Paso 1
    //pause
    //disp("Interaccion numero: "+string(i));
    k1 = funcionLogaritmica(i, L2, ultimoactivo, linealizacion, valMedio, difLimites,nmalla, limites,Coef);
    Laux2(i-1, 2:nmalla+1) = L2(i-1, 2:nmalla+1) + k1(1,1:nmalla) * h/2;
    //pause
    // Paso 2
    k2 =funcionLogaritmica(i, Laux2, ultimoactivo, linealizacion, valMedio, difLimites,nmalla, limites,Coef);
    Laux2(i-1, 2:nmalla+1) = L2(i-1, 2:nmalla+1) + k2(1,1:nmalla) * h/2;
    //pause
    // Paso 3
    k3 = funcionLogaritmica(i, Laux2, ultimoactivo, linealizacion, valMedio, difLimites,nmalla, limites,Coef);
    Laux2(i-1, 2:nmalla+1) = L2(i-1, 2:nmalla+1) + k3(1,1:nmalla) * h;
    //pause
    // Paso 4
    k4 = funcionLogaritmica(i, Laux2, ultimoactivo, linealizacion, valMedio, difLimites,nmalla, limites,Coef);
    L2(i, 2:nmalla+1) = L2(i-1, 2:nmalla+1) + h/6 * (k1(1,1:nmalla) + 2*k2(1,1:nmalla) + 2*k3(1,1:nmalla) + k4(1,1:nmalla));
    diferencial(i+1999,2:nmalla+1)= (k1(1,1:nmalla) + 2*k2(1,1:nmalla) + 2*k3(1,1:nmalla) + k4(1,1:nmalla))/6;
    diferencialLineal(i-1,2:nmalla+1)=(k1(2,1:nmalla) + 2*k2(2,1:nmalla) + 2*k3(2,1:nmalla) + k4(2,1:nmalla))/6;
    diferencialLS(i-1,2:nmalla+1)=(k1(3,1:nmalla) + 2*k2(3,1:nmalla) + 2*k3(3,1:nmalla) + k4(3,1:nmalla))/6;
    diferencialLI(i-1,2:nmalla+1)=(k1(4,1:nmalla) + 2*k2(4,1:nmalla) + 2*k3(4,1:nmalla) + k4(4,1:nmalla))/6;
    diferencialMM(i-1,2:nmalla+1)=(k1(5,1:nmalla) + 2*k2(5,1:nmalla) + 2*k3(5,1:nmalla) + k4(5,1:nmalla))/6;
    diferencialMMLS(i-1,2:nmalla+1)=(k1(6,1:nmalla) + 2*k2(6,1:nmalla) + 2*k3(6,1:nmalla) + k4(6,1:nmalla))/6;
    diferencialMMLI(i-1,2:nmalla+1)=(k1(7,1:nmalla) + 2*k2(7,1:nmalla) + 2*k3(7,1:nmalla) + k4(7,1:nmalla))/6;
    diferencialML(i-1,2:nmalla+1)=(k1(8,1:nmalla) + 2*k2(8,1:nmalla) + 2*k3(8,1:nmalla) + k4(8,1:nmalla))/6;
    diferencialSalida(i-1,2:nmalla+1)=(k1(9,1:nmalla) + 2*k2(9,1:nmalla) + 2*k3(9,1:nmalla) + k4(9,1:nmalla))/6;
    //pause
    total=0;
    // Monitorización
    for m=1:nmalla
        total = total +L2(i, m+1)*valMedio(m);
    end
    totalMasa(i)=total
    if modulo(t, 1000) == 0 then
        mprintf("Iter %d - t=%.2f - Masa=%.2f - y1=%.4f\n", i, t, total, L2(i, 2));
        contIter=contIter+1;
        tiempoActual = getdate("s") - t_inicio; 
        cpu_delta = timer(); // <-- Tiempo CPU de estas últimas 1000 iteraciones
        cpu_acumulado = cpu_acumulado + cpu_delta; // <-- Tiempo CPU total acumulado
        tiempoIter(contIter,1) = i;
        tiempoIter(contIter,2) = tiempoActual;
        tiempoIter(contIter,3) = tiempoActual - tiempoAnt;
        tiempoIter(contIter,4) = cpu_acumulado; 
        tiempoIter(contIter,5) = cpu_delta;
       
        // Actualizamos para la próxima vuelta
        tiempoAnt = tiempoActual;
    end  
    
end
/////////////////////////////////////////////////////////////////////////////
//////////////Distribucion a Malla Uniforme//////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ResultadoTotal=L;
PoblacionTotal(1:n,2:v+1)= L(1:n,2:v+1)*diag(1:v);
PoblacionTotal(1:n,1)=L(1:n,1);
ResultadoTotal(n+1:iterfinal-1,1:linealizacion+1)=L2(2:iterfinal-n,1:linealizacion+1);
PoblacionTotal(n+1:iterfinal-1,2:linealizacion+1)= L2(2:iterfinal-n,2:linealizacion+1)*diag(1:linealizacion);
PoblacionTotal(n+1:iterfinal-1,1)=L2(2:iterfinal-n,1);
q=zeros(iterfinal-n,nmalla);
b=zeros(iterfinal-n,nmalla);
denominador_q2(linealizacion+1:nmalla)=(log10(valMedio(linealizacion+1:nmalla))-log10(valMedio(linealizacion+1-2:nmalla-2)));
denominador_q1(linealizacion+1:nmalla-1)=(log10(valMedio(linealizacion+2:nmalla))-log10(valMedio(linealizacion:nmalla-2)));
denLimites_21(linealizacion+1:nmalla)=limites(linealizacion+1:nmalla,2)-limites(linealizacion:nmalla-1,2);
denLimites_22(linealizacion+1:nmalla)=limites(linealizacion:nmalla-1,1)-limites(linealizacion-1:nmalla-2,1);
denLimites_11(linealizacion+1:nmalla-1)=limites(linealizacion+2:nmalla,2)-limites(linealizacion+1:nmalla-1,2);
denLimites_12(linealizacion+1:nmalla-1)=limites(linealizacion+1:nmalla-1,1)-limites(linealizacion:nmalla-2,1);
cont=0;
for l=linealizacion+1:nmalla-1
    Ceros(1:iterfinal-n)=sign(L2(1:iterfinal-n, l+1))';
    Noceros= find(Ceros==1,1)-1;
    SINO(1:iterfinal-n) = sign(L2(1:iterfinal-n, l+2))';
    limiteInf_1= find(SINO==1,1);
    limiteSup_1= find(SINO==0,1)-1;
    if isempty(limiteSup_1)& limiteInf_1==1 then
        limiteSup_1=length(SINO);
        q(1:iterfinal-n,l)=(log10(L2(1:iterfinal-n,l+2)/denLimites_11(l))-log10(L2(1:iterfinal-n,l)/denLimites_12(l)))/denominador_q1(l);            
    elseif limiteInf_1==1 then
        limiteInf_0=limiteSup_1+1;
        limiteSup_0=length(SINO);
        q(1:limiteSup_1,l)=(log10(L2(1:limiteSup_1,l+2)/denLimites_11(l))-log10(L2(1:limiteSup_1,l)/denLimites_12(l)))/denominador_q1(l);
        q(limiteInf_0:limiteSup_0,l)=(log10(L2(limiteInf_0:limiteSup_0,l+1)/denLimites_21(l))-log10(L2(limiteInf_0:limiteSup_0,l-1)/denLimites_22(l)))/denominador_q2(l);
    elseif isempty(limiteInf_1) & limiteSup_1==0 then
        limiteSup_0=length(SINO);
        q(limiteInf_0:limiteSup_0,l)=(log10(L2(limiteInf_0:limiteSup_0,l+1)/denLimites_21(l))-log10(L2(limiteInf_0:limiteSup_0,l-1)/denLimites_22(l)))/denominador_q2(l);
    elseif limiteInf_1<>1 then
        limiteInf_0= find(SINO==0,1);
        limiteSup_0= find(SINO==1,1)-1;
        limiteSup_1=length(SINO);
        q(limiteInf_0:limiteSup_0,l)=(log10(L2(limiteInf_0:limiteSup_0,l+1)/denLimites_21(l))-log10(L2(limiteInf_0:limiteSup_0,l-1)/denLimites_22(l)))/denominador_q2(l);
        q(limiteInf_1:limiteSup_1,l)=(log10(L2(limiteInf_1:limiteSup_1,l+2)/denLimites_11(l))-log10(L2(limiteInf_1:limiteSup_1,l)/denLimites_12(l)))/denominador_q1(l);
    end
    b(1:iterfinal-n,l)=(q(1:iterfinal-n,l)+1).*L2(1:iterfinal-n,l+1).*(limites(l,2).^(q(1:iterfinal-n,l)+1))./((limites(l,2)+1).^(q(1:iterfinal-n,l)+1)-limites(l,1).^(q(1:iterfinal-n,l)+1));
    for m=limites(l,1):limites(l,2)
                 valInf=m;
                 valSup=m+1;
                 ResultadoTotal(n+1:iterfinal,m+1)=b(1:iterfinal-n,l).*limites(l,1).*((valSup/limites(l,1)).^(q(1:iterfinal-n,l)+1)-(valInf/limites(l,1)).^(q(1:iterfinal-n,l)+1))./(q(1:iterfinal-n,l)+1);
                 PoblacionTotal(n+1:iterfinal,m+1)=ResultadoTotal(n+1:iterfinal,m)*m;
    end
    if Noceros<>0 then
        ResultadoTotal(n+1:n+Noceros,limites(l,1)+1:limites(l,2)+1)=0;
        PoblacionTotal(n+1:n+Noceros,limites(l,1)+1:limites(l,2)+1)=0;
    elseif isempty(Noceros) then
        ResultadoTotal(n+1:iterfinal,limites(l,1)+1:limites(l,2)+1)=0;
        PoblacionTotal(n+1:iterfinal,limites(l,1)+1:limites(l,2)+1)=0;
    end
end
Ceros(1:iterfinal-n)=sign(L2(1:iterfinal-n, nmalla+1))';
Noceros= find(Ceros==1,1)-1;
if isempty(Noceros) then
    ResultadoTotal(n+1:iterfinal,limites( nmalla,1)+1:limites( nmalla,2)+1)=0;
    PoblacionTotal(n+1:iterfinal,limites( nmalla,1)+1:limites( nmalla,2)+1)=0;
else
    q(1:iterfinal-n,nmalla)=q(1:iterfinal-n,nmalla-1);
    b(1:iterfinal-n,nmalla)=(q(1:iterfinal-n,nmalla)+1).*L2(1:iterfinal-n,nmalla+1).*(limites(nmalla,2).^(q(1:iterfinal-n,nmalla)+1))./((limites(nmalla,2)+1).^(q(1:iterfinal-n,nmalla)+1)-limites(nmalla,1).^(q(1:iterfinal-n,nmalla)+1));
    for m=limites(nmalla,1):limites(nmalla,2)
        valInf=m;
        valSup=m+1;
        ResultadoTotal(n+1:iterfinal,m+1)=b(1:iterfinal-n,nmalla).*limites(nmalla,1).*((valSup/limites(nmalla,1)).^(q(1:iterfinal-n,nmalla)+1)-(valInf/limites(nmalla,1)).^(q(1:iterfinal-n,nmalla)+1))./(q(1:iterfinal-n,nmalla)+1);
        PoblacionTotal(n+1:iterfinal,m+1)=ResultadoTotal(n+1:iterfinal,m+1)*m;
    end
    if Noceros<>0 then
        ResultadoTotal(n+1:n+Noceros,limites( nmalla,1)+1:limites( nmalla,2)+1)=0;
        PoblacionTotal(n+1:n+Noceros,limites( nmalla,1)+1:limites( nmalla,2)+1)=0;
    end
end


tiempo= getdate("s") - t_inicio; 
disp("La simulacion ha terminado, en "+string(tiempo))

/////////////////////////////////////////////////////////////////////////////
////////////////////////////Guardar//////////////////////////////////
////////////////////////////////////////////////////////////////////////////
nombreCarpeta = "Log Soule de "+string(no)+" partículas iniciales y de "+string(iterfinal)+" iteraciones don Difusion de Rango z=" +string(z);
ruta_base = "/home/bianca/Documentos/2026/Marzo/";
ruta_carpeta = ruta_base + nombreCarpeta;
mkdir(ruta_carpeta);
csvWrite(PoblacionTotal,ruta_carpeta+"/ResultadoPoblacion.csv")
csvWrite(ResultadoTotal,ruta_carpeta+"/ResultadoConcentracion.csv")
csvWrite(diferencial,ruta_carpeta+"/ResultadoDiferenciales.csv")
csvWrite(limites,ruta_carpeta+"/ResultadoLimites.csv")
csvWrite(linealizacion,ruta_carpeta+"/Linealizacion.csv")
csvWrite(n,ruta_carpeta+"/iteracionesLineales.csv")
csvWrite(valMedio,ruta_carpeta+"/ResultadoValoresMedios.csv")
csvWrite(difLimites,ruta_carpeta+"/ResultadoDifLimites.csv")
csvWrite(diferencialLineal,ruta_carpeta+"/DiferencialLineal.csv")
csvWrite(diferencialLI,ruta_carpeta+"/DiferencialLI.csv")
csvWrite(diferencialLS,ruta_carpeta+"/DiferencialLS.csv")
csvWrite(diferencialML,ruta_carpeta+"/DiferencialML.csv")
csvWrite(diferencialMM,ruta_carpeta+"/DiferencialMM.csv")
csvWrite(diferencialMMLI,ruta_carpeta+"/DiferencialMMLI.csv")
csvWrite(diferencialMMLS,ruta_carpeta+"/DiferencialMMLS.csv")
csvWrite(diferencialSalida,ruta_carpeta+"/DiferencialSalida.csv")
csvWrite(tiempoIter, ruta_carpeta+"/TiempoIteracion.csv")
disp("Termino")
