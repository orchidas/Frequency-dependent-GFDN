function PP = svdPerBin(TT)
    PP = zeros(size(TT,[1 3]));
    for f = 1:size(TT,3)
       s = (svd(TT(:,:,f)));  
       PP(:,f) = s;
    end
end