function [coeffs,rSquare,residualSumSquare,rootMeanSquare] = fitCoeff(cl_a, cl_b, cl_c, dp, dac)
    clustersAcc = (sum(cl_a(dp,:,dac)) + sum(cl_b(dp,:,dac)) + sum(cl_c(dp,:,dac)));
    %=============State=============
    %terminos distintos de 0 en estado
    nonZerosIndexS_a = transpose(find(cl_a(dp,:,dac)~=0));
    nonZerosIndexS_b = transpose(find(cl_b(dp,:,dac)~=0));
    nonZerosIndexS_c = transpose(find(cl_c(dp,:,dac)~=0));
    %largo de clusters asociados a cada estado
    clustersLength = vertcat(nonZerosIndexS_a,nonZerosIndexS_b,nonZerosIndexS_c);
    %frecuencias de clusters
    clustersFreq = vertcat(transpose(cl_a(dp,nonZerosIndexS_a,dac)),transpose(cl_b(dp,nonZerosIndexS_b,dac)),transpose(cl_c(dp,nonZerosIndexS_c,dac)))/clustersAcc;
    %limpiar variables
    clear nonZ* clustersAcc
    
    if length(clustersFreq) >= 2
        y_exp = log(clustersFreq);
        [fitResult, gof] = fit(clustersLength, y_exp, 'poly1');
        coeffs = coeffvalues(fitResult);
        % Coeficiente de determinación R^2
        rSquare = gof.rsquare;
        % Suma de cuadrados de los residuos
        residualSumSquare = gof.sse;
        % Raíz del error cuadrático medio
        rootMeanSquare = gof.rmse;
        
    else
        coeffs = NaN;rSquare = NaN;residualSumSquare = NaN;rootMeanSquare = NaN;
    end
end