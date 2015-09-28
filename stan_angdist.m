function D=stan_angdist(M1,M2)

[ntrials1]=size(M1,2);
[ntrials2]=size(M2,2);

D=zeros(ntrials1,ntrials2);

for i=1:ntrials1
	for j=1:ntrials2
		D(i,j)=median(abs(M1(:,i)-M2(:,j)));
	end
end
