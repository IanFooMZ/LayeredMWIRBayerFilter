classdef functionsContainer
   methods (Static)
       
function fn = formFileName(file,paramCell)
    fc = functionsContainer;
    file{end+1} = fc.compileParamStr(paramCell);

    fn = [];
    for i = 1:length(file)
        if i <= 2% || i == length(file)
            fn = [fn,file{i}];
        else
            fn = [fn,'_',file{i}];
        end
    end
end

function paramStr = compileParamStr(paramCell)
    paramStr = [];
    for k = 1:length(paramCell)/3
        paramStr = [paramStr, '_', paramCell{3*(k-1)+1}, num2str(paramCell{3*(k-1)+2}, paramCell{3*(k-1)+3})];
    end 
    paramStr = paramStr(2:end);
end

function dataArray = extractParamVals(devCell, fieldName)
    dataArray = zeros(length(devCell),length(devCell{1}.sp));
    for k = 1:length(devCell)
        dataArray(k,:) = devCell{k}.(fieldName);
    end
end

function [contrast_processed] = calcContrast(arr_main,arr_1,arr_2,arr_3,peakInd)
    [a, ind] = max(arr_main(peakInd,:));
    arr_processed_main = arr_main(:,ind);
    arr_processed_1 = arr_1(:,ind);
    arr_processed_2 = arr_2(:,ind);
    arr_processed_3 = arr_3(:,ind);
    
    compressArr = [arr_processed_1 arr_processed_2 arr_processed_3];
    nextHighestInt = max(compressArr,[],2);
    
    contrast_processed = arr_processed_main./nextHighestInt;
end

function [Emag_processed] = findMax(Emag,peakInd)
    [a, ind] = max(Emag(peakInd,:));
    Emag_processed = Emag(:,ind);
end

function [power_lambda] = integrateOverSpace(monitor_name,variable,theta)
    Emag = magnitudeField(monitor_name,variable);

    if variable=='E'
        Emag_lambda = squeeze(sum(Emag,[1,2,3]));
        intensity_lambda = 0.5*3e8*8.85418782e-12*1*cos(theta)*Emag_lambda.^2;
        power_lambda = intensity_lambda*(9e-5*9e-5);
    else
        power_lambda = abs(Emag);
    end
end

function [mag_xyz_lambda] = magnitudeField(monitor_name,variable)
    monitor_name = process_dataset(monitor_name,variable);

    field = monitor_name.(variable);
    if variable=='E' || variable=='H'
        fieldX = field(:,:,:,:,1);
        fieldY = field(:,:,:,:,2);
        fieldZ = field(:,:,:,:,3);
        mag_xyz_lambda = sqrt(abs(fieldX).^2+abs(fieldY).^2+abs(fieldZ).^2);
    else
        mag_xyz_lambda = field;
    end
end

function [E_spatial] = reshape_spatial(monitor_name,variable,field_component,wlVal)
%     Deprecated because even though this is probably right I'd rather not
%     deal with this: Ian Foo, 20210924
%     component: x,y,z = 1,2,3
%     wlVal = index in the wavelength value array
    s1 = length(monitor_name.x);
    s2 = length(monitor_name.y);
    s3 = length(monitor_name.z);

    E = monitor_name.(variable)(:,field_component,wlVal);

    E_spatial = reshape(E,[s1 s2 s3]);
end

function [offsetVector] = findOffset(peakInd,valueVectors)
    ind = zeros(length(valueVectors),1);

    for cnt = 1:length(valueVectors)
        [a, ind(cnt)] = max(valueVectors{cnt});
        ind(cnt) = peakInd - ind(cnt);
    end
    
    offsetVector = ind;
end

   end
end
