function [ q,Q,t,outlet ] = read_Q( rivFlx1_file_and_path, riv_file_and_path )

if ( nargin~=2 )
    fprintf('[q(m,1),Q(m,n),time]=readQ(rivFlx1_file_and_path, riv_file_and_path)\n\n\n');
    return
end

flow = load(rivFlx1_file_and_path);

riv=read_riv(riv_file_and_path);
outlet=find(riv(:,4)<0);

t=flow(:,1);
Q=flow(:,2:end);
q=flow(:,outlet+1);

end

