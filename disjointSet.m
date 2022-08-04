classdef disjointSet < handle
    %disjointSet 并查集
    %   find
    %   union

    properties
        set
%         best
    end

    methods
        function obj = disjointSet(best,num)
            % 构造
            %   from best 数组
            obj.set = zeros(1,num);
            for i=1:length(best)
%                 if best(i)
                obj.set(best(i))=best(i);
            end
%             obj.set = best;
%             obj.best = best;
        end

        function index = findDJ(obj,x)
            %findDJ_Naive
            %   递归查询
%             if x>length(obj.set)
%                 x
%             end
%             disp(x)
            if obj.set(x)==x || obj.set(x)==0
                index = obj.set(x);
                return
            end
            obj.set(x) = findDJ(obj,obj.set(x));
            index = obj.set(x);
            return
        end

        function merge(obj, in_1, in_2, neib)
%             merge_count=merge_o+1;
            if obj.set(in_1)*obj.set(in_2)==0
                return
            end
            if length(neib{in_1})<length(neib{in_2})
                [in_1, in_2] = deal(in_2,in_1); % in_1.neib >= in_2.neib
            end
            obj.set(findDJ(obj,in_2))=findDJ(obj,in_1);
            obj.set(in_1) = findDJ(obj,obj.set(in_1));
            obj.set(in_2) = findDJ(obj,obj.set(in_2));
        end
    end
end