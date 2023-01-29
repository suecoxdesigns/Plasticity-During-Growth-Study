function d2 = extractJump2019c(sublist,j,d2)
%recalculate the weight
weights = [sublist.weight];
weights = weights(weights>2);
averageWeight = nanmean(weights);
d2(j).weight =averageWeight;
d2(j).dates = unique({sublist.date2});
d2(j).AniNum = sublist(1).AniNum;
if (d2(j).AniNum == 762 || d2(j).AniNum == 788 || d2(j).AniNum==800)
    d2(j).AniNum
end

ones = find([sublist.jumpMaxNorm]==1);
if (~isempty(ones) &&(~isempty(averageWeight)&&~isnan(averageWeight)))
    [sublist(ones).jumpMaxZero] = deal(averageWeight);
end

d2(j).Treat = sublist(1).Treat;
d2(j).AniNumS =  sublist(1).AniNumS;
if d2(j).AniNum == 761
    j=j;
end

% combine all the jumps into one struct
jumps = [sublist.jump];

% for each jump - calculate stuff
if ~isempty(jumps)
    for i=1:length(jumps)
        jumps(i).impulseNonZeroed = sum(jumps(i).jumpFz);
        % calcuate power
        dt = jumps(i).jumpT(2)-jumps(i).jumpT(1);
        %first calculate vertical acceleration
        jumps(i).accel = jumps(i).jumpFz./(averageWeight/9.8);
        jumps(i).v(1) = 0;
        for k=2:length(jumps(i).accel)
            jumps(i).v(k) = jumps(i).v(k-1)+dt*(jumps(i).accel(k)+jumps(i).accel(k-1))/2;
        end
        jumps(i).power = jumps(i).jumpFz'.*jumps(i).v;
        jumps(i).maxPower = max([jumps(i).power]);
        jumps(i).maxV = max([jumps(i).v])
        
    end
    d2(j).maxImpulseNonZeroed = max([jumps(i).impulseNonZeroed]);
    d2(j).maxPower = max([jumps.maxPower]);
    d2(j).maxV = max([jumps.maxV]);
    
else
    d2(j).maxImpulseNonZeroed =0;
    d2(j).maxPower = 0;
    d2(j).maxV;
end
d2(j).AllJumps = [sublist.jumpMaxZero];
d2(j).AllJumps = d2(j).AllJumps(d2(j).AllJumps>1);
trialTimeSecTotal = sum([sublist.trialTimeSec]);
d2(j).AllImpulseZ = [sublist.zImpulses];
d2(j).nJumps = length(d2(j).AllJumps);
d2(j).jumpPsec = d2(j).nJumps/trialTimeSecTotal;
%subset out just the jumps and calculate the corefficient of variation
jumpsAllInfo = [sublist.jump];

if ~isempty(jumpsAllInfo)
    %%
    [~,order] = sort([jumpsAllInfo.jumpMaxNorm],'descend');
    if length(order)>2
        jumpsAllInfoTop3 = jumpsAllInfo(order(1:3));
    else
        jumpsAllInfoTop3 = jumpsAllInfo(order);
    end
    clear idxStart idxEnd top3fzNorm meanTop3Fz sdTop3Fz power
    top3fzNorm = zeros(length(jumpsAllInfoTop3),(101));
   
    for i=1:length(jumpsAllInfoTop3)
        [~,idxStart] = min(abs([jumpsAllInfoTop3(i).jumpTLong]-(-0.1)));
        [val,idxEnd] = min(abs([jumpsAllInfoTop3(i).jumpTLong]-(0.07)));
        tempFz = jumpsAllInfoTop3(i).jumpFzLong(idxStart:idxEnd)';
        dataIn=tempFz;
        top3fzNorm(1,:) = (spline(1:1:length(dataIn),dataIn,0:length(dataIn)/100:(length(dataIn))));
        jumptime(i) = jumpsAllInfoTop3(i).jumpT(end)-jumpsAllInfoTop3(i).jumpT(1);
     %  power(i) = jumpsAllInfoTop3(i).impulse^2/(2*averageWeight/9.8)/jumptime(i);
    end
    
    meanTop3Fz = nanmean(top3fzNorm,1);
    sdTop3Fz = std(top3fzNorm,0,1);
    if length(jumpsAllInfoTop3)>1
        d2(j).sdAv = (sum(sdTop3Fz.^2)/length(sdTop3Fz))^(0.5);
        d2(j).CoefVarAv = 100*length(sdTop3Fz)*d2(j).sdAv/(sum(abs(meanTop3Fz)));
    end
    d2(j).meanJumpTimeTop3 = nanmean([jumptime]);
    [d2(j).maxImpulse,idx] = nanmax([jumpsAllInfoTop3.impulse]);
    d2(j).AvTop3Impulses = nanmean([jumpsAllInfoTop3.impulse]);
    d2(j).jumpTimeTopImpulse = jumptime(idx);
   % d2(j).maxPower = nanmax(power);
    %recalculate the normalized Jumps
    if ~isempty(averageWeight)
        d2(j).AllJumpsNorm = d2(j).AllJumps/averageWeight;
        
    else
        d2(j).AllJumpsNorm = [sublist.JumpMaxNorm];
    end
    JumpMaxZTot = sum([d2(j).AllJumpsNorm]);
    d2(j).MaxJumpTot = nanmax([d2(j).AllJumpsNorm]);
    d2(j).forceTimeTot = JumpMaxZTot/trialTimeSecTotal;
    
    % average of top 3 jumps
    clear jumps
    [jumps,order] = sort([d2(j).AllJumpsNorm],'descend');
    

    if length(jumps)>2
        d2(j).Avoftop3MaxJumps = nanmean(jumps(1:3));
        d2(j).TopJumps = jumps(1:3);
    elseif length(jumps)==2
        d2(j).Avoftop3MaxJumps = nanmean([jumps,1]);
        d2(j).TopJumps = jumps(1:2);
    elseif length(jumps)==1
        d2(j).Avoftop3MaxJumps = nanmean([jumps,1,1]);
        d2(j).TopJumps = jumps;
    end
    if isempty(d2(j).MaxJumpTot)
        d2(j).Avoftop3MaxJumps = 1;
        d2(j).forceTimeTot = 0;
        d2(j).MaxJumpTot= 1;
        d2(j).TopJumps = 1;
    end
end
end