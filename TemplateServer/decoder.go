package main

// Adapted from https://raw.githubusercontent.com/dogenzaka/tsv/master/parser.go
// MIT License

import (
	"encoding/csv"
	"encoding/json"
	"errors"
	"io"
	"reflect"
	"strconv"

	"github.com/DisposaBoy/JsonConfigReader"
	"github.com/go-playground/validator/v10"
)

// Parser has information for parser
type Parser struct {
	Headers    []string
	Reader     *csv.Reader
	Data       interface{}
	ref        reflect.Value
	indices    []int // indices is field index list of header array
	structMode bool
}

// NewTsvParser creates new TSV parser with given io.Reader
func NewTsvParser(reader io.Reader, data interface{}) *Parser {
	r := csv.NewReader(reader)
	r.Comma = '\t'

	p := &Parser{
		Reader: r,
		Data:   data,
		ref:    reflect.ValueOf(data).Elem(),
	}

	return p
}

// Next puts reader forward by a line
func (p *Parser) Next() (eof bool, err error) {
	// Get next record
	var records []string

	for {
		// read until valid record
		records, err = p.Reader.Read()
		if err != nil {
			if err.Error() == "EOF" {
				return true, nil
			}
			return false, err
		}
		if len(records) > 0 {
			break
		}
	}

	if len(p.indices) == 0 {
		p.indices = make([]int, len(records))
		// mapping simple index
		for i := 0; i < len(records); i++ {
			p.indices[i] = i + 1
		}
	}

	// record should be a pointer
	for i, record := range records {
		idx := p.indices[i]
		if idx == 0 {
			// skip empty index
			continue
		}
		// get target field
		field := p.ref.Field(idx - 1)
		switch field.Kind() {
		case reflect.String:
			field.SetString(record)
		case reflect.Bool:
			if record == "" {
				field.SetBool(false)
			} else {
				col, err := strconv.ParseBool(record)
				if err != nil {
					return false, err
				}
				field.SetBool(col)
			}
		case reflect.Float64, reflect.Float32:
			if record == "" {
				field.SetFloat(0)
			} else {
				col, err := strconv.ParseFloat(record, 64)
				if err != nil {
					return false, err
				}
				field.SetFloat(col)
			}
		case reflect.Int64, reflect.Int32, reflect.Int16, reflect.Int8, reflect.Int:
			if record == "" {
				field.SetInt(0)
			} else {
				col, err := strconv.ParseInt(record, 10, 64)
				if err != nil {
					return false, err
				}
				field.SetInt(col)
			}
		case reflect.Uint64, reflect.Uint32, reflect.Uint16, reflect.Uint8, reflect.Uint:
			if record == "" {
				field.SetUint(0)
			} else {
				col, err := strconv.ParseUint(record, 10, 64)
				if err != nil {
					return false, err
				}
				field.SetUint(col)
			}
		default:
			return false, errors.New("Unsupported field type")
		}
	}

	return false, nil
}

func DecodeJsonAndValidate(r io.Reader, target interface{}) error {
	dec := json.NewDecoder(JsonConfigReader.New(r))

	err := dec.Decode(target)
	if err != nil {
		return err
	}

	validate := validator.New()
	if err := validate.Struct(target); err != nil {
		return err
	}

	return nil
}
