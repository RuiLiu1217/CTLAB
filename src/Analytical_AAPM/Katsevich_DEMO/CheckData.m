function Res = CheckData(Data,DataMin, DataMax)

Res=Data;
if(Data<DataMin)
    Res= DataMin;
end;
if(Data>DataMax)
    Res= DataMax;
end;
