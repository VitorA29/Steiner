void UseFPmax(DadosMD& dadosMD){
	int n_padroes = 10;
	ostringstream buffer;
	buffer.str("");		
//                   fpmax_hnmp <semente> <id_arq_tmp> <banco de dados> <tam. do banco> <suporte minimo> <qtd de padroes> <arq. saida>

	buffer << "./fpmax_hnmp " << "1 " << (genrand_int32() % 100) << " bd.txt " <<  dadosMD.listaSolucoesElite.size() << " " << dadosMD.freq_supp << " " << n_padroes << " padroes.txt" ;
//	buffer << "./fim_maximal " << "bd.txt " << dadosMD.freq_supp << " padroes.txt" ;
	int s = system(buffer.str().c_str());
	buffer.str("");
		
}