function MakeTitle(fid, titulo, larguraTotal)
    numPontos = (larguraTotal - length(titulo)) / 2;
    linhaPontos = repmat(':', 1, floor(numPontos));
    linhaTitulo = [linhaPontos titulo linhaPontos];
    fprintf(fid,'\n');
    fprintf(fid,linhaTitulo);
    fprintf(fid,'\n');
end