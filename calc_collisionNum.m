function numInter = calc_collisionNum(Floe)

numInter=cat(1,Floe.interactions);
if ~isempty(numInter)
    
    numInter=size(numInter(numInter(:,1)<Inf,1),1)/2+size(numInter(numInter(:,1)==Inf,1),1);
    
else
    numInter=0;
end

end

