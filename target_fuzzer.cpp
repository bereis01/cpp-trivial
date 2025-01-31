#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <string>

using namespace std;

#define NUM_AMINOS 20
#define PENALIDADE_INDEL 0

vector<vector<int>> matriz_pontuacao_blosum62 = {   { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
                                                    {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
                                                    {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
                                                    {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
                                                    { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
                                                    {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
                                                    {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
                                                    { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
                                                    {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
                                                    {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
                                                    {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
                                                    {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
                                                    {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
                                                    {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
                                                    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
                                                    { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
                                                    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
                                                    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
                                                    {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
                                                    { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};


unordered_map<char, int> dicionario_indice_alfabeto_all_amino = {
    {'A', 0},
    {'R', 1},
    {'N', 2},
    {'D', 3},
    {'C', 4},
    {'Q', 5},
    {'E', 6},
    {'G', 7},
    {'H', 8},
    {'I', 9},
    {'L', 10},
    {'K', 11},
    {'M', 12},
    {'F', 13},
    {'P', 14},
    {'S', 15},
    {'T', 16},
    {'W', 17},
    {'Y', 18},
    {'V', 19}
};

pair<vector<vector<int>>, vector<vector<int>>> needleman_wunsch_iterativo(
    const string& v,
    const string& w,
    const vector<vector<int>>& matriz_pontuacao = matriz_pontuacao_blosum62,
    const unordered_map<char, int>& dicionario_indice_alfabeto = dicionario_indice_alfabeto_all_amino,
    int penalidade_indel = PENALIDADE_INDEL)
{
    int len_v = v.size();
    int len_w = w.size();

    // A matriz que armazena as regras de posicionamento ou direcionamento do caminho
    // 0 = "north", 1 = "west" e 2 = "northwest"
    vector<vector<int>> matriz_posicionamento(len_v + 1, vector<int>(len_w + 1, 0));
    fill(matriz_posicionamento[0].begin(), matriz_posicionamento[0].end(), 1);

    // Inicializando a matriz de recorrência ou memoização
    vector<vector<int>> matriz_similaridade(len_v + 1, vector<int>(len_w + 1, 0));
    for (int i = 1; i <= len_v; ++i) {
        matriz_similaridade[i][0] = matriz_similaridade[i - 1][0] - penalidade_indel;
    }

    for (int j = 1; j <= len_w; ++j) {
        matriz_similaridade[0][j] = matriz_similaridade[0][j - 1] - penalidade_indel;
    }

    for (int i = 1; i <= len_v; ++i) {
        for (int j = 1; j <= len_w; ++j) {
            // Inserção em v:
            int valor_insercao_v = matriz_similaridade[i - 1][j] - penalidade_indel;

            // Inserção em w:
            int valor_insercao_w = matriz_similaridade[i][j - 1] - penalidade_indel;

            // Casamento
            int valor_match = matriz_similaridade[i - 1][j - 1] + matriz_pontuacao[dicionario_indice_alfabeto.at(v[i - 1])][dicionario_indice_alfabeto.at(w[j - 1])];

            // Atualizando o valor na matriz de similaridade
            matriz_similaridade[i][j] = max({valor_insercao_v, valor_insercao_w, valor_match});

            // Atualizando a operação que deve ser feita no caminhamento
            if (matriz_similaridade[i][j] == valor_match) {
                matriz_posicionamento[i][j] = 2;
            } else if (matriz_similaridade[i][j] == valor_insercao_w) {
                matriz_posicionamento[i][j] = 1;
            } else if (matriz_similaridade[i][j] == valor_insercao_v) {
                matriz_posicionamento[i][j] = 0;
            }
        }
    }

    return {matriz_posicionamento, matriz_similaridade};
}

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *Data, size_t Size) {
    try
    {
        // Transforma a entrada do Fuzzer em uma string
        string data_as_str(reinterpret_cast<const char*>(Data), Size);

        size_t meio = Size / 2;

        // Quebra a string no meio. v = primeira metade e w = segunda metade
        needleman_wunsch_iterativo(data_as_str.substr(0, meio), data_as_str.substr(meio));
    }
    catch(const std::exception& e)
    {
        return 0;
    }
    
	return 0;
}