%


% nbird, ndays, nrois

morning=[];
evening=[];
between=[];

frac=4;

for i=1:length(stats)
  for j=1:length(stats(i).vmat)

    ntrials=size(stats(i).vmat{j},2);
    pool1=1:round(ntrials/frac);
    pool2=round((ntrials-ntrials/frac)):ntrials;

    for k=1:size(stats(i).vmat{j},3)
      morning=[morning mean(mean(stats(i).vmat{j}(pool1,pool1,k)))];
      evening=[evening mean(mean(stats(i).vmat{j}(pool2,pool2,k)))];
      between=[between mean(mean(stats(i).vmat{j}(pool1,pool2,k)))];
    end
  end
end
