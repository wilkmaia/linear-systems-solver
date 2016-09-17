extern crate time;

use std::io::prelude::*;
use std::fs::File;


/// Ler do arquivo de entrada
/// Lê a entrada formatada do arquivo 'input.txt'
/// buf: String auxiliar para a leitura do arquivo
/// output: Vetor de strings onde cada elemento do vetor é uma linha do arquivo de entrada
fn read_input_file<'a>( mut buf: &'a mut String ) -> Vec<&'a str> {
	// Tentativa de leitura do arquivo de entrada
	let mut f = match File::open( "input.txt" ) {
		Ok( file ) => file,
		Err( _ ) => panic!( "Erro ao abrir o arquivo de entrada." ),
	};
	
	// Armazenamento do conteúdo do arquivo no buffer
	match f.read_to_string( &mut buf ) {
		Ok( _ ) => (), // Se a leitura tiver sido bem sucedida, o resultado estará em "buf"
		Err( _ ) => panic!( "Erro ao ler o arquivo de entrada." ),
	};
	
	// Retorna com um vetor organizado na forma
	// ["Linha 1", "Linha 2", "Linha 3", ...]
	buf.split( "\r\n" ).collect::<Vec<&str>>()
}


/// Escrever no Arquivo de Saída
/// Escreve a saída formatada como pedido no arquivo 'output.txt'
/// x: Vetor solução do sistema
/// r: Resíduos do sistema
/// dif: Precisão obtida na saída do sistema
/// z: Quantidade de iterações necessárias para métodos iterativos
/// n: Dimensão do vetor solução
/// method: Método utilizado
/// p: Precisão numérica utilizada
fn write_to_output_file( x: Vec<f64>, r: Vec<f64>, dif: Vec<f64>, z: u32, n: usize, method: u8, p: usize ) {
	// Tentativa de abertura do arquivo de saída
	let mut f = match File::create( "output.txt" ) {
		Ok( file ) => file,
		Err( _ ) => panic!( "Erro ao criar o arquivo de saída." ),
	};
	
	let mut s = String::new();
	let mut i: usize = 0;

	// Formata o vetor solução
	while i < n {
		s.push_str( &( format!( "{:.*}", p, x[i] ) ) );
		s.push_str( "\r\n" );
		i += 1;
	}
	
	// Formata o vetor de resíduos
	i = 0;
	while i < n {
		s.push_str( &( format!( "{:.*}", p, r[i] ) ) );
		s.push_str( "\r\n" );
		i += 1;
	}
	
	// Formata o vetor da precisão obtida
	i = 0;
	while i < n {
		s.push_str( &( format!( "{:.*}", p, dif[i] ) ) );
		s.push_str( "\r\n" );
		i += 1;
	}

	// Caso seja um método iterativo, exibe a quantidade de iterações
	if method == 2  || method == 3 {
		s.push_str( &( z.to_string() ) );
	}
	
	// Tentativa de escrita no arquivo de saída
	match f.write_all( s.as_bytes() ) {
		Ok( _ ) => (),
		Err( _ ) => panic!( "Erro ao escrever no arquivo de saída." ),
	};
}


/// Ajuste de Matrizes
/// Ajusta a matriz A de tal forma que os maiores elementos de cada
/// coluna estejam na diagonal principal, evitando elementos nulos e
/// números em ponto flutuante muito grandes
/// a: Matriz de coeficientes
/// b: Vetor de resultados
/// n: Dimesão da matriz A
fn ajuste_de_matrizes<'a>( mut a: &'a mut Vec<Vec<f64>>, mut b: &'a mut Vec<f64> ) {
	// Dimensão da matriz x
	let n = b.len();

	// Varredura da matriz a fim de ordenar as linhas para evitar
	// problemas nos cálculos. Técnica de Pivoteamento Parcial.
	let mut i: usize = 0;
	while i < n {
		// Encontra o maior valor da coluna, em módulo
		let mut _max: f64 = a[i][i].abs();
		let mut lin: usize = i;
		let mut j: usize = i + 1;
		while j < n {
			if _max < a[j][i].abs() {
				lin = j;
				_max = a[j][i].abs();
			}

			j += 1;
		}

		// Caso haja necessidade de efetuar a troca entre linhas
		// Isso acontece aqui
		if lin != i {
			let _aux: Vec<f64> = a[i].to_vec();
			a[i]	= a[lin].to_vec();
			a[lin]	= _aux.to_vec();

			_max	= b[i];
			b[i]	= b[lin];
			b[lin]	= _max;
		}

		i += 1;
	}
	// Fim da remoção dos elementos nulos da diagonal principal
}


/// Método de Gauss
/// Aplica a técnica de Pivoteamento ou Pivoteamento Parcial para a resolução
/// de sistemas lineares do tipo ax = b
/// Args
/// a: Matriz de coeficientes
/// b: Vetor de resultados
/// x: Incógnitas do problema
fn metodo_de_gauss<'a>( mut a: &'a mut Vec<Vec<f64>>, mut b: &'a mut Vec<f64>, mut x: &'a mut Vec<f64> ) {
	// Dimensão da matriz x
	let n = x.len();

	let mut k: usize = 0;
	let mut i: usize = 1;
	let mut z: usize = 0;
	let mut _k: usize = 0;
	
	// Quantidade de manipulações necessárias
	let count = ( n - 1 ) * n / 2;

	while z < count {
		k = if i < n { k } else { i = k + n - 2; k + 1 };

		// Aplicação do método
		let pivot = a[k][k];
		let m = - a[i][k] / pivot;
		
		// Atualização da linha 'i'
		let mut _i = 0;
		while _i < n {
			a[i][_i] += m * a[k][_i];
			_i += 1;
		}
		b[i] += m * b[k];
		
		z += 1;
		i += 1;
	}
	
	// Com a matriz a triangularizada, calculam-se os valores
	// de x[i] de forma direta, com i = (n-1)..0
	i = n;
	while i > 0 {
		i -= 1;

		let mut sum: f64 = 0f64;
		k = n;
		while k > i {
			k -= 1;

			sum += x[k] * a[i][k];
		}
		x[i] = ( b[i] - sum ) / a[i][i];
	}
}


/// Métodos Iterativos
/// Aplica os métodos de Jacobi ou de Gauss-Seidel para a resolução de sistemas lineares
/// Args
/// a: Matriz de coeficientes
/// b: Vetor de resultados
/// x: Incógnitas do problema
/// dif: Vetor com a precisão dos resultados alcançados
/// z: Quantidade de iterações necessárias para a convergência do método
/// er: Precisão exigida
/// method: Método a ser aplicado (2 -> Jacobi / 3 -> Gauss-Seidel)
fn metodos_iterativos<'a>( a: &'a Vec<Vec<f64>>, b: &'a Vec<f64>, mut x: &'a mut Vec<f64>, mut dif: &'a mut Vec<f64>, mut z: &'a mut u32, er: f64, method: u8 ) {
	// Dimensão da matriz x
	let n = x.len();

	// Chute inicial -> 0
	let mut y = vec![0f64; n];
	*z = 0;
	
	// Execução indefinida do laço até a convergência
	loop {
		let mut i: usize = 0;

		// Cálculo da resposta para a i-ésima iteração
		while i < n {
			let mut sum: f64 = 0f64;
			let mut j = 0;

			// Soma das parcelas a[i][j] * x[j]
			while j < n {
				// O termo a[i][i] será o denominador da expressão
				// Não entra na soma
				if i == j {
					j += 1;
					continue;
				}

				// A única diferença entre os métodos de Jacobi e Gauss-Seidel é que
				// a cada iteração, os valores atualizados são utilizados, quando disponíveis.
				//
				// Isto quer dizer que para a etapa i1 > i0 desta soma, o valor de x a ser utilizado
				// é o valor calculado na mesma iteração (y[j]), ao invés de x[j], que representa 
				// a iteração anterior.
				sum = sum + a[i][j] * ( if method == 2 { x[j] } else { y[j] } );

				j += 1;
			}
			
			// Atualização do valor de X(i) da iteração atual
			y[i] = ( b[i] - sum ) / a[i][i];

			i += 1;
		}
		
		i = 0;
		let mut flag = true;

		// Cálculo da precisão alcançada na iteração atual
		while i < n {
			dif[i] = ( x[i] - y[i] ).abs();
			if dif[i] > er {
				flag = false;
			}

			i += 1;
		}

		// Atualização do vetor X
		*x = y.to_vec();

		*z += 1;

		// Caso o sistema tenha convergido
		if flag == true {
			break;
		}
	}
}


/// Função principal
/// Coordena a entrada do sistema, etapas de computação e saída
fn main() {
	// Leitura do arquivo de entrada
	let mut buf = String::new();
	let lines: Vec<&str> = read_input_file( &mut buf );
	
	// Parsing dos dados
	// Método a ser utilizado
	let method = lines[0].parse::<u8>().unwrap();
	
	// Dimensão da matriz A
	let n = lines[1].parse::<usize>().unwrap();
	
	// Precisão exigida para os métodos iterativos
	let er = lines[2].parse::<f64>().unwrap();
	
	// Matriz A
	let mut i: usize = 0;
	let mut a: Vec<Vec<f64>> = Vec::with_capacity( n*n );
	while i < n {
		let numbers = lines[3+i].split( " " );
		let numbers = numbers.collect::<Vec<&str>>();
		let mut j = 0;
		let mut aux: Vec<f64> = Vec::with_capacity( n );
		while j < n {
			aux.push( match numbers[j].parse::<f64>() {
				Ok( num ) => num,
				Err( _ ) => panic!( "Problema na leitura da matriz A." ),
			} );
			j += 1;
		}
		a.push( aux );
		i += 1;
	}
	
	// Matriz B
	let mut k: usize = 0;
	let mut b: Vec<f64> = Vec::with_capacity( n );
	while k < n {
		b.push( match lines[3+i+k].parse::<f64>() {
			Ok( num ) => num,
			Err( _ ) => panic!( "Problema na leitura da matriz B." ),
		} );
		k += 1;
	}
	// Fim do parsing dos dados

	
	// Matriz X
	let mut x = vec![0f64; n];
	
	// Auxiliares
	let mut z: u32 = 0u32;
	let mut dif: Vec<f64> = vec![0f64; n];
	
	// Execução do método proposto
	// O tempo de execução é avaliado em nanossegundos
	let t1 = time::precise_time_ns();

	ajuste_de_matrizes( &mut a, &mut b );

	match method {
		1 => metodo_de_gauss( &mut a, &mut b, &mut x ),
		2|3 => metodos_iterativos( &a, &b, &mut x, &mut dif, &mut z, er, method ),
		_ => panic!( "Método inválido!" ),
	};

	let dt = time::precise_time_ns() - t1;
	// Fim da execução do método proposto e da medição do tempo gasto
	
	// Exibição do tempo de execução e do método aplicado
	println!( "Tempo de execução: {} us", dt/1000 );
	println!( "Método: {}", method );

	
	// Cálculo do vetor de resíduos
	let mut r: Vec<f64> = Vec::with_capacity( n );
	let mut sum: f64;
	i = 0;
	while i < n {
		let mut j = 0;
		sum = 0f64;
		while j < n {
			sum += a[i][j] * x[j];
			j += 1;
		}
		r.push( b[i] - sum );
		i += 1;
	}

	// println!( "{:?}", a ); // DEBUG
	// println!( "{:?}", b ); // DEBUG


	// Escrita no arquivo de saída
	write_to_output_file( x, r, dif, z, n, method, ( ( 10f64 / er ).log10() as usize ) );
	// ( 10f64 / er ).log10() as usize
	// Utilizado para passar a precisão dos resultados numéricos a serem exibidos
	// no arquivo de saída.
	//
	// A quantidade de casas decimais exibidas é 1 + N, onde
	// N = log10( precisão^(-1) )
	//
	// Logo, para uma precisão de 0.001, serão exibidas
	// N = log10( 1000 ) = 3
	// (1 + 3) = 4 casas decimais

	// write_to_output_file( x, r, dif, z, n, method, 50 as usize ); // DEBUG
}
