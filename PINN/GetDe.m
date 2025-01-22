function De = GetDe(K, G)

De = zeros(6, 6);
De(1,1) = K+4/3*G;
De(2,2) = K+4/3*G;
De(3,3) = K+4/3*G;
De(4,4) = G;
De(5,5) = G;
De(6,6) = G;
De(1,2) = K-2/3*G;
De(1,3) = K-2/3*G;
De(2,1) = K-2/3*G;
De(2,3) = K-2/3*G;
De(3,1) = K-2/3*G;
De(3,2) = K-2/3*G;

end