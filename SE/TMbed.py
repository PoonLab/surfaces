import sys
import torch
import argparse
import os
from tqdm import tqdm
from pathlib import Path
from tmbed.model import Predictor
from tmbed.embed import T5Encoder
from tmbed.viterbi import Decoder
from tmbed.utils import seed_all, read_fasta, make_batches, collate_batch, make_mask

#sys.path.append('/content/tmbed/')
#seed_all(101)

models = load_models(Path('/home/sareh/tmbed/tmbed/models/cnn/'))
encoder = load_encoder(Path('/home/sareh/tmbed/tmbed/models/t5/'))
decoder = Decoder()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str,
                       help='dir to aa fasta files')
    parser.add_argument('--outdir', type=str,
                       help='path to directory to write outfile')                   
    return parser.parse_args()


def load_encoder(model_path):
    return T5Encoder(model_path, use_gpu=False)


def load_models(model_path):
    models = []

    for model_file in sorted(model_path.glob('*.pt')):
        model = Predictor()

        model.load_state_dict(torch.load(model_file)['model'])

        model = model.eval()

        models.append(model)

    return models


def predict_sequences(models, embeddings, mask):
    B, N, _ = embeddings.shape

    num_models = len(models)

    with torch.no_grad():
        pred = torch.zeros((B, 5, N), device=embeddings.device)

        for model in models:
            y = model(embeddings, mask)
            pred = pred + torch.softmax(y, dim=1)

        pred = pred / num_models

    return pred.detach()


def write_3_line(output_file, proteins, predictions, pred_map):
    with output_file.open('w') as of:
        for protein in proteins:
            seq_hash = protein.seq_hash

            if seq_hash not in predictions:
                continue

            header = protein.header
            sequence = protein.sequence

            prediction, probabilities = predictions[seq_hash]

            prediction = ''.join(pred_map[v] for v in prediction.tolist())

            assert len(prediction) == len(sequence)

            of.write(f'{header}\n')
            of.write(f'{sequence}\n')
            of.write(f'{prediction}\n')


def write_tabular(output_file, proteins, predictions, pred_map):
    col_head = 'AA\tPRD\tP(B)\tP(H)\tP(S)\tP(i)\tP(o)'

    with output_file.open('w') as of:
        for protein in proteins:
            seq_hash = protein.seq_hash

            if seq_hash not in predictions:
                continue

            header = protein.header
            sequence = protein.sequence

            prediction, probabilities = predictions[seq_hash]

            prediction = ''.join(pred_map[v] for v in prediction.tolist())

            probabilities = probabilities.tolist()

            assert len(prediction) == len(sequence)
            assert len(probabilities) == len(prediction)

            of.write(f'{header}\n')
            of.write(f'{col_head}\n')

            for aa, prd, probs in zip(sequence, prediction, probabilities):
                probs = '\t'.join(f'{v:.2f}' for v in probs)

                of.write(f'{aa}\t{prd}\t{probs}\n')

def main():
    args = parse_args()
    for file in os.listdir(args.indir):
        infile = os.path.join(args.indir, file)
        outfile = os.path.join(args.outdir, file)
    
        proteins = read_fasta(infile)
        encoder.to_cuda()
        device = encoder.device()

        for i in range(len(models)):
            models[i] = models[i].to(device)

        predictions = dict()

        sorted_proteins = sorted(proteins, key=lambda protein: protein.length)

        batches = make_batches(sorted_proteins, batch_size=4000)

        with tqdm(total=len(sorted_proteins), leave=True) as progress:
            for a, b in batches:
                batch = sorted_proteins[a:b]

                lengths = [protein.length for protein in batch]
                sequences = [protein.sequence for protein in batch]

                embeddings = encoder.embed(sequences)

                embeddings = embeddings.to(device)
                embeddings = embeddings.to(dtype=torch.float32)

                mask = make_mask(embeddings, lengths)

                probabilities = predict_sequences(models, embeddings, mask)

                mask = mask.cpu()
                probabilities = probabilities.cpu()

                prediction = decoder(probabilities, mask).byte()

                probabilities = probabilities.permute(0, 2, 1)

                for idx, protein in enumerate(batch):
                    length = protein.length
                    seq_hash = protein.seq_hash
                    predictions[seq_hash] = (prediction[idx, :length],
                                             probabilities[idx, :length])

                progress.update(b - a)

        encoder.to_cpu()

        for i in range(len(models)):
            models[i] = models[i].cpu()

        torch.cuda.empty_cache()

        # --out-format 0 or 2
        pred_map = {0: 'B', 1: 'b', 2: 'H', 3: 'h', 4: 'S', 5: '.', 6: '.'}

        # --out-format 0 or 1
        write_3_line(output_file, proteins, predictions, pred_map)

if __name__ == "__main__":
    main()
