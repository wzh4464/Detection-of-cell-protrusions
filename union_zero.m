function new_ssets=union_zero(ssets,neib)
    if size(ssets,2)<=1
        new_ssets = ssets;
        return
    end
    for i=1:size(ssets,2)
        for j=2:size(ssets,2)
            for k=1:size(ssets{i})
                for m=1:size(ssets{j})
                    if find(neib{ssets{j}(m)}==ssets{i}(k))
                        temp_1=horzcat(ssets{i},ssets{j});
                        ssets{i}=unique(temp_1);
                        ssets(j)=[];
                        new_ssets=union_zero(ssets,neib);
                        return
                    end
                end
            end
        end
    end
    new_ssets=ssets;
end