function A= fastpowermaxplus(B, t)
    exp = bin(t)
    value = B
    for i =3:lenght(exp)
        value= otimes(value,value)
        print i,":\t",value
        if(exp(i:i+1)=='1')
            value = otimes(value,B)
	   
        end
    end
end