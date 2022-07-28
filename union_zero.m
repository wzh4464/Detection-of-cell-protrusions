function new_ssets=union_zero(ssets,neib)
    if size(ssets,2)<=1
        new_ssets = ssets;
        return
    end
    for i=1:size(ssets,2) % ssets{i}
        for j=i+1:size(ssets,2) % ssets{j}
            for k=1:size(ssets{i})
                for m=1:size(ssets{j})
                    if find(neib{ssets{j}(m)}==ssets{i}(k)) %ssets{j}.neib intersect ssets{i} is not null
                        temp_1=horzcat(ssets{i},ssets{j});
                        ssets{i}=unique(temp_1);
                        ssets(j)=[];
                        disp([size(ssets,2)])
                        new_ssets=union_zero(ssets,neib);                    
                        return
                    end
                end
            end
        end
    end
    new_ssets=ssets;
end