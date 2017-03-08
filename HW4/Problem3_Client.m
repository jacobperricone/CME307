clear all
close all

i
npoints = 1000
x = generate_x_points(npoints)
s = generate_s_points(npoints)


p = [6,12]
pname = {'six', 'twelve'}
cutoffs = [0,-1];
cname = {'zero', 'negTen'};



Values = zeros(4, npoints, 4)


k = 1
for i=1:size(p,2)
    for j=1:size(cutoffs,2)
        [out, point] = evaluate_candidacy(x,s,p(i), cutoffs(j));
        Values(1,:,k) = out
        Values(2:4,1:size(point,2),k) = point
        k = k + 1 
    end
    
  
    
end

%% p = 6

p6 = Values(:,:, 1:2)
outputs6 = squeeze(p6(1,:,:))
points6 = p6(2:4,:,:)

figure()
subplot(2,1,1)


%%
p12 = Values(:,:, 3:4)
outputs12 = squeeze(p12(1,:,:));
points12 = p12(2:4,:,:)



