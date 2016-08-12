%Makes arrays of coordinates from "origin" to final", according to
%"directions"

fid = fopen( 'stations.txt', 'wt' );

dim = 3;
freq = 2;
Vs = 5/sqrt(2);
lambda = freq*Vs;


origin = [5.0,5.0,5.0];
final  = [50.0,50.0,50.0];
diff = final - origin;
directions = [1, 0, 0; 
             0, 1, 0; 
             0, 0, 1; 
             1, 1, 1];
         
nLambda = max(diff)/lambda;
nPoints = nLambda*10;

vec = linspace(0,1,nPoints+1)'; %+1 For the origin
Mat = repmat(vec,1,dim);


%Writing Center
fprintf(fid, '%15.6f ', origin);
fprintf(fid, '\n');

%Writing others
for i = 1:size(directions,1)
    array = Mat;
    %dire  = directions(i,:);
    dire  = directions(i,:)/norm(directions(i,:));
    for d = 1:dim
        array(:,d) = (dire(d)*array(:,d)*diff(d))+origin(d);
    end
    
    for n=2:size(array,1)
       fprintf(fid, '%15.6f ', array(n,:));
       fprintf(fid, '\n');
    end
  
end

    

% Mat = zeros(numel(vec), dim);
% Mat = Mat+origin;
% 
% %Axis
% for i = 1:dim
%     Mat = zeros(numel(vec), dim);
%     Mat = Mat+origin;
%     
%     Mat(:,i) = Mat(:,i) + vec;
%     a=2;
% end
% 
% %Diagonal
% Mat = zeros(numel(vec), dim);
% Mat = Mat+origin;
%  for i = 1:dim   
%     Mat(:,i) = Mat(:,i) + vec;
%     a=2;
% end   
    
%for image = 1:N
%  [a1,a2,a3,a4] = ProcessMyImage( image );
%  fprintf( fid, '%f,%f,%f,%f\n', a1, a2, a3, a4);

fclose(fid);