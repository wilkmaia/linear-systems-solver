#![allow(non_snake_case)]

extern crate time;

use std::io::prelude::*;
use std::fs::File;

fn read_input_file<'a>( mut buf: &'a mut String ) -> Vec<&'a str> {
	let mut f = match File::open( "input.txt" ) {
		Ok( file ) => file,
		Err( _ ) => panic!( "Erro ao abrir o arquivo de entrada." ),
	};
	
	match f.read_to_string( &mut buf ) {
		Ok( _ ) => (), // Se a leitura tiver sido bem sucedida, o resultado estará em "buf"
		Err( _ ) => panic!( "Erro ao ler o arquivo de entrada." ),
	};
	
	// Retorna com um vetor organizado na forma
	// ["Linha 1", "Linha 2", "Linha 3", ...]
	buf.split( "\r\n" ).collect::<Vec<&str>>()
}


fn write_to_output_file( x: Vec<f64>, r: Vec<f64>, dif: Vec<f64>, z: u32, n: usize, method: u8, p: usize ) {
	let mut f = match File::create( "output.txt" ) {
		Ok( file ) => file,
		Err( _ ) => panic!( "Erro ao criar o arquivo de saída." ),
	};
	
	let mut s = String::new();
	let mut i: usize = 0;
	while i < n {
		s.push_str( &( format!( "{:.*}", p, x[i] ) ) );
		s.push_str( "\r\n" );
		i += 1;
	}
	
	i = 0;
	while i < n {
		s.push_str( &( format!( "{:.*}", p, r[i] ) ) );
		s.push_str( "\r\n" );
		i += 1;
	}
	
	i = 0;
	while i < n {
		s.push_str( &( format!( "{:.*}", p, dif[i] ) ) );
		s.push_str( "\r\n" );
		i += 1;
	}
	if method == 2  || method == 3 {
		s.push_str( &( z.to_string() ) );
	}
	
	match f.write_all( s.as_bytes() ) {
		Ok( _ ) => (),
		Err( _ ) => panic!( "Erro ao escrever no arquivo de saída." ),
	};
}


/// Método de Gauss
/// Aplica a técnica de Pivoteamento ou Pivoteamento Parcial para a resolução
/// de sistemas lineares do tipo ax = b
/// Args
/// a: Matriz de coeficientes
/// b: Vetor de resultados
/// x: Incógnitas do problema
/// pivoteamento_parcial: epecifica se deve ser utilizada a técnica de
///						  pivoteamento (false) ou pivoteamento_parcial (true)
fn metodo_de_gauss<'a>( mut a: &'a mut Vec<Vec<f64>>, mut b: &'a mut Vec<f64>, mut x: &'a mut Vec<f64>, pivoteamento_parcial: bool ) {
	// Dimensão da matriz x
	let n = x.len();

	let mut k: usize = 0;
	let mut i: usize = 1;
	let mut z: usize = 0;
	let mut _k: usize = 0;
	
	// Quantidade de manipulações necessárias
	let count = ( n - 1 ) * n / 2;

	// Método de Pivoteamento utilizado
	while z < count {
		k = if i < n { k } else { i = k + n - 2; k + 1 };

		// Ajuste das linhas - Pivoteamento Parcial
		// Realizado apenas uma vez para cada valor de k (linha base)
		if k == _k && pivoteamento_parcial == true {
			let mut j = k + 1;
			let mut _max: f64 = a[k][k].abs();
			let mut lin: usize = k;
			while j < n { // Verifica qual a linha com maior pivô
				if _max < a[j][k].abs() {
					lin = j;
					_max = a[j][k].abs();
				}

				j += 1;
			}

			// Verifica se há necessidade de troca de linhas
			// Caso haja, troca-as
			if lin != k {
				let mut _aux: Vec<f64> = Vec::with_capacity( n );

				_aux	= a[k].to_vec();
				a[k]	= a[lin].to_vec();
				a[lin]	= _aux.to_vec();

				_max	= b[lin];
				b[lin]	= b[k];
				b[k]	= _max;
			}

			_k = k + 1;
			println!("Ajuste de linhas: k = {} // i = {}", k, i);
		}
		// Fim do ajuste de linhas - Pivoteamento Parcial

		// Aplicação do método
		let pivot = a[k][k];
		let m = - a[i][k] / pivot;
		
		let mut _i = 0;
		while _i < n {
			a[i][_i] = m * a[k][_i] + a[i][_i];
			_i += 1;
		}
		
		b[i] = m * b[k] + b[i];
		
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


fn metodos_iterativos<'a>( a: &'a Vec<Vec<f64>>, b: &'a Vec<f64>, mut x: &'a mut Vec<f64>, mut dif: &'a mut Vec<f64>, mut z: &'a mut u32, er: f64, method: u8 ) {
	// Dimensão da matriz x
	let n = x.len();

	// Chute inicial -> 0
	let mut y = vec![0f64; n];
	*z = 0;
	
	loop {
		let mut i: usize = 0;

		while i < n {
			let mut sum: f64 = 0f64;
			let mut j = 0;

			while j < n {
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
				sum += a[i][j] * ( if method == 2 { x[j] } else { y[j] } );

				j += 1;
			}
			
			y[i] = ( b[i] - sum ) / a[i][i];
			i += 1;
		}
		
		i = 0;
		let mut flag = true;

		while i < n {
			dif[i] = ( x[i] - y[i] ).abs();
			if dif[i] > er {
				flag = false;
			}

			i += 1;
		}

		// println!("x = {:?}", x);
		// println!("y = {:?}", y);
		// println!("");
		
		*x = y.to_vec();
		*z += 1;
		
		// Caso o sistema tenha convergido
		if flag == true {
			break;
		}
	}
}


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
	
	// Matriz X
	let mut x = vec![0f64; n];
	
	// Auxiliares
	let mut z: u32 = 0u32;
	let mut dif: Vec<f64> = vec![0f64; n];
	
	// Execução do método proposto
	// O tempo de execução é avaliado em nanossegundos
	let t1 = time::precise_time_ns();
	
	match method {
		1 => metodo_de_gauss( &mut a, &mut b, &mut x, true ),
		2|3 => metodos_iterativos( &a, &b, &mut x, &mut dif, &mut z, er, method ),
		_ => panic!( "Método inválido!" ),
	};

	let dt = time::precise_time_ns() - t1;
	
	println!( "Tempo de execução: {} ns", dt );
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
}
