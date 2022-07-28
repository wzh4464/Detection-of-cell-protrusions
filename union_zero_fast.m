function [connectSet,num]=union_zero_fast(best,neib)
%% neib to edge
% %     edge=zeros(length(neib),2);
%     edges = cell(1,length(neib)); 
%     parfor i=1:length(neib)
% %         edge(i)=[];
%         for j=1:length(neib{i})
%             if find(best==i)
%                 break
%             end
%             if i<neib{i}(j)
%                 edges{i} = [edges{i};neib{i}(j)];
%             end
%         end
%     end % edges{i} 表示(i,edges{i})线段
%     merge_count=0;
    connectDS=disjointSet(best);
    for i=1:length(best) % merge best中每个值和其edge邻接
        for j=1:length(neib{best(i)})
%             disp([i,j])
            connectDS.merge(best(i),neib{best(i)}(j),neib);
        end
    end
%     for i=1:length(best) % merge best中每个值和其edge邻接
%         con
%     end
    parfor i=1:length(neib)
        stan(i)=i
    end
    num=size(find(connectDS.set==stan),2);
    connectSet=connectDS.set;
%     edge = unique(edge);
%%
%     connectSet=[];
%     for i=1:length(neib)
%         connectSet=[connectSet;edge{i}];
%     end
end