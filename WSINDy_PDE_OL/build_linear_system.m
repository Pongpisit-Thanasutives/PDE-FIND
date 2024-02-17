function [G,b,M]=build_linear_system(Theta_cell,Cfs_t_scaled,lib_list,lhs_ind,M_full)

dim = length(size(Theta_cell{1}));
foo=repmat({1},dim-1,1);
lib_list_cell=mat2cell(lib_list,ones(size(lib_list,1),1),size(lib_list,2));
Theta_pdx = cell2mat(cellfun(@(x,y)reshape( convn(x,reshape(Cfs_t_scaled(y(end)+1,:),foo{:},[]),'valid') ,[],1),...
    Theta_cell,lib_list_cell,'uni',0)');
num_eq = length(lhs_ind);
[K,m] = size(Theta_pdx);
G = Theta_pdx(:,~ismember(1:m,lhs_ind));
b = zeros(K,num_eq);
M = [];
for k=1:num_eq
    b(:,k) = Theta_pdx(:,lhs_ind(k));
    if ~isempty(M_full)
        M = [M M_full(~ismember(1:m,lhs_ind))/M_full(lhs_ind(k))];
    end
end