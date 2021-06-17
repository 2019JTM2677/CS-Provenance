function index = FindIndex(n,a,b,dest_id)

    if b == dest_id
        index = (n-1)*a;
    elseif a < b
        index = (n-1)*(b-1)+a;
    else 
        index = (n-1)*(b-1)+a-1;
    end
end