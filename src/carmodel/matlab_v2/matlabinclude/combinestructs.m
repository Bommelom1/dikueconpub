function p1 = combinestructs(varargin)
    % COMBINESTRUCTS: Procedure that combines of several structurs struct1, struct1, struct3,.... 
    %
    % SYNTAX:
    %   p1 = struct2vec(struct1, struct1, struct3,...)
    %
    % OUTPUTS
    %   pstruct:    Structure that combines information form input structurs struct1, struct1, struct3,.... 
    %               
    %               If input structures have fields with the same names, the fields in the first arguments will
    %               be recusively overwritten/ubdated with values in the later arguments.   
    % 
    % INPUTS
    %   struct1, struct2.. :  Structures with fields taking scalar values (at least one field for each value element in pnames) 
    % 
    % Date November 4
    
    p1=varargin{1};
    for j=2:nargin; 
        p2=varargin{j};
        pnames=fieldnames(p2);       
        k=numel(pnames);
        for i=1:k;
            p1.(char(pnames(i)))=p2.(char(pnames(i)));
        end
    end
end % end of combinestructs