package main

import (
	"bufio"
	"io"
	"math"
	"os"
	"sort"
)

type Entry struct {
	Key    string
	Offset uint64
	Length uint64
}

type Index []Entry

func (a Index) Len() int           { return len(a) }
func (a Index) Swap(i, j int)      { a[i], a[j] = a[j], a[i] }
func (a Index) Less(i, j int) bool { return a[i].Key < a[j].Key }

type Reader struct {
	Index Index
	file  *os.File
}

func (d *Reader) Make(data string, index string) {
	file, err := os.Open(data)
	if err != nil {
		return
	}
	d.file = file

	f, err := os.Open(index)
	if err != nil {
		return
	}
	defer f.Close()

	d.Index = make([]Entry, 0)
	entry := Entry{}
	parser := NewTsvParser(bufio.NewReader(f), &entry)
	for {
		eof, err := parser.Next()
		if eof {
			break
		}
		if err != nil {
			return
		}
		d.Index = append(d.Index, entry)
	}
	sort.Sort(Index(d.Index))
}

func (d *Reader) Delete() {
	d.file.Close()
}

func (d *Reader) Id(key string) (int64, bool) {
	i := int64(sort.Search(int(d.Size()), func(i int) bool {
		return d.Index[i].Key >= key
	}))
	return i, i < d.Size() && d.Index[i].Key == key
}

func (d *Reader) Key(id int64) (string, bool) {
	if id < 0 || id >= d.Size() {
		return "", false
	}
	return d.Index[id].Key, true
}

func (d *Reader) Offset(id int64) uint64 {
	if id < 0 || id >= d.Size() {
		return math.MaxUint64
	}
	return d.Index[id].Offset
}

func (d *Reader) Length(id int64) uint64 {
	if id < 0 || id >= d.Size() {
		return math.MaxUint64
	}
	return d.Index[id].Length
}

func (d *Reader) Data(id int64) string {
	if id < 0 || id >= d.Size() {
		return ""
	}
	d.file.Seek(int64(d.Index[id].Offset), io.SeekStart)
	length := d.Index[id].Length - 1
	buffer := make([]byte, length)
	d.file.Read(buffer)
	return string(buffer[:length])
}

func (d *Reader) Size() int64 {
	return int64(len(d.Index))
}
