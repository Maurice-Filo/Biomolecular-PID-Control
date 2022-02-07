function Struct = CombineStructures(Struct1, Struct2)
Struct = Struct1;
FN_Struct2 = fieldnames(Struct2);
for i = 1:length(FN_Struct2)
	Struct.(FN_Struct2{i}) = Struct2.(FN_Struct2{i});
end  
end

